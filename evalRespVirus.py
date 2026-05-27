"""
Script for regular evaluation and overview of respiratory viral infections.

:author: Christopher Hardt <christopher.hardt@laborberlin.com>

:sources:
- Leistungsverzeichnis: http://s-labb-mds01.laborberlin.com/csmed/App/MDS/Catalog
- https://www.arcgis.com/apps/dashboards/2a29f2bebc524c67b6250b64beea12bf
- https://public.data.rki.de/t/public/views/ARE-Dashboard/BevoelkerungGrippeWeb-Inzidenzen?%3Aembed=y&%3AisGuestRedirectFromVizportal=y

"""

import numpy as np
import pandas as pd
import os
import sys
sys.path.insert(0, '..')
sys.path.insert(0, '/home/chardt/projects/BID/bid-korlab/')
import korlab as kl
import matplotlib.pyplot as plt
from matplotlib.dates import MonthLocator, DateFormatter
import pyexasol
from dotenv import load_dotenv
from datetime import date, timedelta
import json


import smtplib #for sending the email
# Import the email modules we'll need
#from email.message import EmailMessage
from email.mime.text import MIMEText
from email.mime.image import MIMEImage
from email.mime.multipart import MIMEMultipart


targetResult = 'positive'
validResults = {
    'neg': 'negative',
    'pos': 'positive',
    'nina': 'negative',
    'hpos': 'positive',
    'mpos': 'positive',
    'gpos': 'positive',
    'stpos': 'positive',
    'swpos': 'positive',
}


# PUBLICATION SETTINGS
LIVE = True # whether to send email to VIR and BID; if False: only send to chardt (for testing)
SEND_MAIL = True # whether to send the report as an email with attached plot; if False: only create plot and save on smb
RKI = True # whether to retrieve and integrate current incidence data from the Robert-Koch Institute

# REPORT THRESHOLD SETTINGS
TOPSPECS = ['SARS-CoV-2', 'Influenza A', 'Influenza B', 'RSV'] # Agents that are always reported
MINRATE = 3 # Current rate in % has to be at least this high to be marked as non-green and be reported (default: 3%)
MINPERCENTILE = 'p60' # Current rate has to be above agent's YTD median to be reported

if LIVE:
    EMAIL_FROM = 'bioinformatik@laborberlin.com'
    EMAIL_TO = ['victor.corman@laborberlin.com']
    EMAIL_CC = ['christopher.hardt@laborberlin.com'] #['bioinformatik@laborberlin.com'] 
else:
    EMAIL_FROM = 'no-reply@laborberlin.com'
    EMAIL_TO = ['christopher.hardt@laborberlin.com'] #['christopher.hardt@laborberlin.com'] ['bioinformatik@laborberlin.com']
    EMAIL_CC = []

### RKI sources (incidence data) ###
REPORT_URL = 'https://raw.githubusercontent.com/robert-koch-institut/GrippeWeb_Daten_des_Wochenberichts/refs/heads/main/GrippeWeb_Daten_des_Wochenberichts.tsv'
URLS = {
    'SARS-CoV-2': 'https://raw.githubusercontent.com/robert-koch-institut/COVID-19_7-Tage-Inzidenz_in_Deutschland/refs/heads/main/COVID-19-Faelle_7-Tage-Inzidenz_Bundeslaender.csv',
    'Influenza A': 'https://raw.githubusercontent.com/robert-koch-institut/Influenzafaelle_in_Deutschland/refs/heads/main/IfSG_Influenzafaelle.tsv',
    'RSV': 'https://raw.githubusercontent.com/robert-koch-institut/Respiratorische_Synzytialvirusfaelle_in_Deutschland/refs/heads/main/IfSG_RSVfaelle.tsv'}
DISEASE = 'ARE'
REGION = 'Osten'

OUTDIR = '/home/chardt/smb/KORLAB/VIR'
SPECS = {'cvdp': 'SARS-CoV-2', 'Influenza A': 'Influenza A', 'Influenza B': 'Influenza B', 'rsvp': 'RSV'}
#SPECS.update({'novp': 'Norovirus', 'rtvp': 'Rotavirus', 'mpvp': 'Mumpsvirus', 'msvp': 'Masernvirus', 'vzvp': 'Varizella-Zoster-Virus'})
#SPECS.update({'bpsp': 'Bordetella (para-)pertussis','ruvp': 'Rubellavirus'})
#SPECS.update({'havp': 'Hepatitis-A-Virus', 'hbvp': 'Hepatitis-B-Virus', 'hcvp': 'Hepatitis-C-Virus', 'hevp': 'Hepatitis-E-Virus'})
#SPECS.update({'mytp': 'Tuberkulose', 'legp': 'Legionella sp.'})
#SPECS = {'Influenza A': 'Influenza A', 'Influenza B': 'Influenza B'} #FOR TESTING Influenza-only
#SPECS = {'Influenza B': 'Influenza B'} 

#Mapping of Influenza LOINCs to InfA or InfB
#LOINC2INFTYPE = {'34487-9': 'Influenza A', '40982-1': 'Influenza B'}
#UUNR2INFTYPE = {'3081': 'MPV', '18541': 'MPV'}

MULTIPLEX = False
MULTIPLEX_UUNRS = ['18531','3751']
MULTIPLEX_MAP = {'Adenovir': 'Adenovirus', 'CoroHKU1': 'Coronav. endem.', 'CoroOC43': 'Coronav. endem.', 'CoroNL63': 'Coronav. endem.',
'EntRhino': 'Entero-/Rhinovirus', 'InfA09H1': 'Influenza A', 'InflA': 'Influenza A', 'InflB': 'Influenza B', 'MPV': 'MPV',
'Parainf1': 'Parainfluenza', 'Parainf2': 'Parainfluenza', 'Parainf3': 'Parainfluenza', 'Parainf4': 'Parainfluenza',
'Rhinovir': 'Entero-/Rhinovirus', 'RSV': 'RSV', 'SARS2': 'SARS-CoV-2'}

PERIOD = 'ytd'
PLOT = True
ROLL_DAYS = 7
TICK_DIST = 7
TREND_COLORS = {7: 'red', 30: 'gold'} #verbose: {3: 'red', 7: 'gold', 30: 'darkblue', 100: 'purple'}


#Get credentials
if False: #via dotenv
    load_dotenv()
    EXASOLAR_USERuid_e = os.getenv('EXASOL_UID')
    EXASOLAR_PWpw_e = os.getenv('EXASOL_PW')
else:
    #via KeePassXC -> environment variables
    EXASOLAR_USER = 'BIOINFORMATIK_RO'
    EXASOLAR_PW = os.environ['PWD']

def get_rki_data(startDate: str, endDate: str) -> pd.DataFrame:
    df = pd.read_csv(REPORT_URL, sep='\t')
    #print(df.tail())
    #print(df['Region'].value_counts())
    df = df.loc[(df['Erkrankung'] == DISEASE) & (df['Region'] == REGION)] #filter for ARE and East Germany
    #print(df.head())
    #Init empty data frame with index from starting date to end date
    df = pd.DataFrame(index=pd.date_range(start=startDate, end=endDate, freq="D"))
    for spec, url in URLS.items():
        if spec == 'SARS-CoV-2': #daily format
            df_tmp = pd.read_csv(url, sep=',')
            df_tmp = df_tmp.loc[(df_tmp['Bundesland_id'] == 11) & (df_tmp['Altersgruppe'] == '00+')]
            df_tmp.index = pd.to_datetime(df_tmp['Meldedatum'])
            df_tmp = df_tmp['Inzidenz_7-Tage'] #.rename(spec)
            #print(df_tmp.head)
        elif spec in ['Influenza A','RSV']: #weekly format
            df_tmp = pd.read_csv(url, sep='\t')
            df_tmp = df_tmp.loc[(df_tmp['Region_Id'] == 11) & (df_tmp['Altersgruppe'] == '00+')]
            df_tmp.index = pd.to_datetime(df_tmp['Meldewoche'] + "-7", format="%G-W%V-%u")
            df_tmp = df_tmp['Inzidenz'] #.rename(spec)
            #print(df_tmp.head)
        df[spec] = df_tmp #.reindex(df.index, fill_value=0) #.reindex(df.index)
        #exit()
        #df_rsv = df_rsv.loc[(df_rsv['Datum'] >= pd.to_datetime(startDate)) & (df_rsv['Datum'] <= pd.to_datetime(endDate))]
        #df_rsv.rename(columns={'Inzidenz': 'rsvp'}, inplace=True)

    '''
    #Covid-19
    df_rsv = pd.read_csv(URLS['cvdp'], sep='\t')
    df_rsv = df_rsv.loc[(df_rsv['Region_Id'] == 11) & (df_rsv['Altersgruppe'] == '00+')]
    df_rsv['Datum'] = pd.to_datetime(df_rsv['Meldewoche'] + "-7", format="%G-W%V-%u")
    df_rsv = df_rsv[['Datum', 'Inzidenz']]
    df_rsv = df_rsv.loc[(df_rsv['Datum'] >= pd.to_datetime(startDate)) & (df_rsv['Datum'] <= pd.to_datetime(endDate))]
    df_rsv.rename(columns={'Inzidenz': 'rsvp'}, inplace=True)
    print(df_rsv.head())
    #RSV
    df_rsv = pd.read_csv(URLS['rsvp'], sep='\t')
    df_rsv = df_rsv.loc[(df_rsv['Region_Id'] == 11) & (df_rsv['Altersgruppe'] == '00+')]
    df_rsv['Datum'] = pd.to_datetime(df_rsv['Meldewoche'] + "-7", format="%G-W%V-%u")
    df_rsv = df_rsv[['Datum', 'Inzidenz']]
    df_rsv = df_rsv.loc[(df_rsv['Datum'] >= pd.to_datetime(startDate)) & (df_rsv['Datum'] <= pd.to_datetime(endDate))]
    df_rsv.rename(columns={'Inzidenz': 'rsvp'}, inplace=True)

    df_count = df_count.asfreq(freq='D', fill_value=0) #Convert time series to daily, adding 0 for days without measurement
    '''
    print(df.tail(n=20))
   
    return df

def getData(startDate: str, endDate: str) -> pd.DataFrame:
    targetUUNRs = [item for lst in uunr_map.values() for item in lst]
    #print(targetUUNRs)
    #exit()
    #Connect to DB (before: exasol.laborberlin.intern/nocertcheck:8563)
    C = pyexasol.connect(dsn='exasol.laborberlin.intern/nocertcheck', user=EXASOLAR_USER, password=EXASOLAR_PW, schema='PRODUCTION', compression=True)
    #Get how often a species was found per day
    if MULTIPLEX:
        sql = f'''
            SELECT gmwstr, CAST(anfotim AS DATE) AS datum, uunr FROM messwert
            WHERE uunr IN ('{"', '".join(MULTIPLEX_UUNRS)}')
            AND anfotim > '{startDate}'
            AND anfotim < '{endDate}'
        '''
        df = C.export_to_pandas(sql)
        print(df['GMWSTR'].value_counts()[:30])
        exit()
    else:
        sql = f'''
            SELECT m.gmwstr, CAST(m.anfotim AS DATE) AS datum, u.uunr FROM messwert m
            INNER JOIN analyt u ON m.uunr = u.uunr
            WHERE u.uunr IN ('{"', '".join(targetUUNRs)}')
            AND anfotim > '{startDate}'
            AND anfotim < '{endDate}'
        '''


        sql_t = f'''
            SELECT m.gmwstr, CAST(m.anfotim AS DATE) AS datum, u.uunr, u.erregerprofilcode, u.loinc FROM messwert m
            INNER JOIN analyt u ON m.uunr = u.uunr
            WHERE u.erregerprofilcode IS NOT NULL
            AND anfotim > '{startDate}'
            AND anfotim < '{endDate}'
        '''

        #print(sql)
        df = C.export_to_pandas(sql)
        #print(df['UUNR'].value_counts()[:30])
        print(df['GMWSTR'].value_counts()[:30])
        #print(df['ERREGERPROFILCODE'].value_counts())
        #exit()

        uunr2inftype = {int(item): key for key, items in UUNR_MAP.items() for item in items}
        print(uunr2inftype)
        #print(df.head())
        #df['ERREGER'] = df['LOINC'].map(LOINC2INFTYPE).fillna(df['ERREGERPROFILCODE']) #.replace(LOINC2INFTYPE)
        df['ERREGER'] = df['UUNR'].map(uunr2inftype)
        print(df.head(n=100).to_string())
        print(df['ERREGER'].value_counts())
        #exit()
        df['GMWSTR'] = df['GMWSTR'].map(validResults).fillna('invalid')
        print(df['GMWSTR'].value_counts())
        #exit()
        
    df['DATUM'] = pd.to_datetime(df['DATUM'])
    #print(df.info())
    return df

def trend_sentence(
    s: pd.Series,
    trend_window=30,
    min_points=3,
    smooth_window=7
    
) -> tuple[str, dict()]:
    """
    s: pd.Series of values indexed by datetime-like index
    returns: (sentence, diagnostics_dict)
    """

    if isinstance(s, pd.DataFrame):
        s = s.iloc[:, 0]

    # Smooth first by 7d rolling mean
    s_smooth = s.rolling(smooth_window, min_periods=1).mean()

    # Select recent window
    end = s_smooth.index.max() #commonly today's date
    start = end - pd.Timedelta(days=trend_window) #e.g. 30 days back
    w = s_smooth.loc[s_smooth.index >= start].dropna() #only look at time period

    if len(w) < min_points:
        return "Not enough data", {}
    
    # Regression
    x = ((w.index - w.index[0]).total_seconds() / 86400).astype(float) #seconds per day
    y = w.to_numpy(dtype=float)
    b, a = np.polyfit(x, y, 1) #least squares polynomial of degree 1 (= straight)

    #print(b, a)

    # Median and percentiles
    p20, p40, median_level, p60, p80 = np.nanpercentile(s.to_numpy(dtype=float), [20, 40, 50, 60, 80])
    if np.isnan(p20):
        print(np.isnan(s.values).sum())
        print(np.isfinite(s.values).sum())
        print(s.values.dtype)
        print(s.values)
        print(p20)
        exit()
    METHOD = 'pctchg' #polyslope, pctchg
    if METHOD == 'polyslope':
        if b > 0.1:
            direction = "strongly increasing"
        elif b < -0.1:
            direction = "strongly decreasing"
        elif b > 0.05:
            direction = "increasing"
        elif b < -0.05:
            direction = "decreasing"
        else:
            direction = "flat"
    elif METHOD == 'pctchg':
        # Percent change over window
        
        pct_change = (b * (x.max() - x.min())) / median_level * 100
        #print(b, x.max(), x.min(), median_level, pct_change)
        #print(f'% change: {pct_change}')
        # Direction
        if pct_change > 100:
            direction = "strongly increasing"
        elif pct_change < -100:
            direction = "strongly decreasing"
        elif pct_change > 20:
            direction = "increasing"
        elif pct_change < -20:
            direction = "decreasing"
        else:
            direction = "flat"
    #print(direction)

    # today's rolling mean level vs YTD median and distribution
    #print(s.iloc[-1])
    today = float(s_smooth.iloc[-1])
    
    textColor = 'green'
    
    if (today < MINRATE) or (today <= p20):
        level = "at a low level"
        textColor = 'green'
    elif today >= p80:
        level = "at a high level"
        textColor = 'red'
    elif today <= p40:
        level = "at a decreased level"
        textColor = '#aaaaaa'
    elif today >= p60:
        level = "at an increased level"
        textColor = 'orange'
    else:
        level = "at a typical level"
        textColor = 'black'

    sentence = f'{direction}, {level}'
    print(today, textColor)
    return sentence, {'current_rate_roll': today, 'pct_change_over_window': pct_change, 'a': float(a), 'b': float(b), 's': start, 'e': end, 'color': textColor, 'p40': p40, 'p50': median_level, 'p60': p60}


def main():
    #end = date.today()
    start, end = kl.getPeriod(PERIOD)
    periods = TREND_COLORS.keys() #[3, 7, 30, 100] #in days

    # Get incidence data from RKI
    if RKI:
        df_rki = get_rki_data(start, end)
        print(df_rki.head())
        #exit()
    # Create the container (outer) email message.
    msg = MIMEMultipart()
    header = [str(x)+'d' for x in periods]
    seenUUNRs = {}
    htmlBody = f'''
        <h2>Aktuelle Entwicklung der Positivraten</h2><table><tr><th style="padding: 5px;">Erreger</th><th style="padding: 5px;">{'</th><th style="padding: 5px;">'.join(header)}</th><th style="padding: 5px;">Pos.Rate (Median&#177;)</th></tr><tr>
    '''

   
    #Get data
    with open('/home/chardt/projects/BID/bid-korlab/VIR/erreger_uunrs.json', 'r', encoding="utf-8") as f:
        uunr_map = json.load(f)

    #uunr_map = {k: uunr_map[k] for k in ['Coronav. endem.']} #FOR TESTING
    uunr_map = {k: uunr_map[k] for k in uunr_map.keys() if ' all' not in k} #filter out specs ending with "all"
    targetUUNRs = [item for lst in uunr_map.values() for item in lst]
    specs = list(uunr_map.keys())
    print(specs)

    uunr2inftype = {int(item): key for key, items in uunr_map.items() for item in items}
    df = kl.getRespVir(EXASOLAR_PW, start, end, targetUUNRs)
    #print(df[['UUNR', 'GMWSTR']].value_counts().to_string())

    df['ERREGER'] = df['UUNR'].map(uunr2inftype)
    #print(df['GMWSTR'].value_counts())

    #Start plot
    fig, axs = plt.subplots(len(specs), 1, sharex=False, figsize=(24, 4*len(specs)))
    plt.xlabel('Datum')
    plotCounter = 0
    
    
    #make 5 groups
    for spec in specs:
        if spec != 'MPV':
            pass #continue
        print(f'Adding {spec}...')
        #fileAtts = f'{spec}_{PERIOD}'
        #df = getData(spec, start, end)
        df_spec = df[df['ERREGER'] == spec]
        #print(df_spec.head())
        #print(df_spec['GMWSTR'].value_counts())
        seenUUNRs[spec] = df_spec['UUNR'].unique().astype('str').tolist()
        
        df_spec = df_spec.drop(['UUNR', 'ERREGER'], axis=1)

        df_pos = df_spec[df_spec['GMWSTR'] == targetResult].groupby('DATUM').count().sort_index()
        df_pos = df_pos.asfreq(freq='D', fill_value=0) #Convert time series to daily, adding 0 for days without measurement
        print(df_pos.tail(n=30))
        #df_sum = df_spec.groupby('DATUM').count().sort_index()
        #df_pos = df_spec[df_spec['GMWSTR'] == targetResult]
        #df_count = df_pos.groupby('DATUM').count().sort_index()
        #df_count = df_count.asfreq(freq='D', fill_value=0) #Convert time series to daily, adding 0 for days without measurement


        #Plot rolling means of absolute counts     
        count_smooth = df_pos.rolling(7, min_periods=1).mean()
        axs[plotCounter].plot(df_pos.index, count_smooth, label="Absolut (7d-Mittel)", linewidth=2, color='deepskyblue')      
        axs[plotCounter].set_ylabel('Positivmessungen (abs.)', color='deepskyblue')
        axs[plotCounter].set_title(spec)
        axs[plotCounter].grid(False)
        axs[plotCounter].set_xlim(np.datetime64(start), np.datetime64(end))
        axs[plotCounter].xaxis.set_major_locator(MonthLocator())
        axs[plotCounter].xaxis.set_major_formatter(DateFormatter('%Y-%m'))

        #Plot rolling means of relative counts (positive rates)
        nontargetResult = 'negative'

        df_neg = df_spec[df_spec['GMWSTR'] == nontargetResult].groupby('DATUM').count().sort_index()
        df_neg = df_neg.asfreq(freq='D', fill_value=0) #Convert time series to daily, adding 0 for days without measurement

        df_total = df_pos.add(df_neg, fill_value=0)
        #df_rate = df_pos / df_total * 100
        df_rate = df_pos.div(df_total, fill_value=0).mul(100) #avoid division by zero, keep NaN if no measurements at all
        rate_smooth = df_rate.rolling(7, min_periods=2).mean()
        '''
        print(rate_smooth.tail(n=30))
        if spec in TOPSPECS:
            if rate_smooth['GMWSTR'].iloc[-1] < MINRATE:
                axs[plotCounter].set_title(spec + f' (aktuell exkludiert, da Positivrate < {MINRATE}%)')
            elif MINPERCENTILE == 'median' and rate_smooth['GMWSTR'].iloc[-1] <= rate_smooth['GMWSTR'].median():
                axs[plotCounter].set_title(spec + f' (aktuell exkludiert, da Positivrate <= Jahresmedian von {rate_smooth["GMWSTR"].median():.1f}%)')
        '''
        ax2 = axs[plotCounter].twinx()
        ax2.plot(df_rate.index, rate_smooth, label="Positivrate (7d-Mittel)", linewidth=2, color='darkorange')      
        #ax2.set_ylabel('Positivrate (in %)', color='darkorange', rotation=270)
        ax2.grid(True)
        ax2.text(
                1.02, 0.5, 'Positivrate (in %)',
                transform=ax2.transAxes,
                rotation=270,
                va='center',
                ha='left',
                color='darkorange'
            )
        

        if not MULTIPLEX: 
            specs = [targetResult]

        # TRENDS (as line plots as well as text)
        htmlTrends = ''
        col = 'black'
        med = 0
        p60 = 0
        for p in periods:
            sentence, params = trend_sentence(df_rate, p)
            print(f'##### {spec} ({p}d)#####\n{params}\n{sentence}')
            if params:
                col = params['color']
                arrow = '&#8680' #rightwards arrow ("flat")
                if 'decreasing' in sentence:
                    arrow = '&#8681' #downwards arrow
                elif 'increasing' in sentence:
                    arrow = '&#8679' #upwards arrow
                if 'strongly' in sentence:
                    arrow += arrow # Add second arrow for strong trends
                htmlTrends += f'<td style="color:{col}; padding: 5px;">{arrow}</td>'
                curRate = params['current_rate_roll']
                
                p40 = params['p40']
                p60 = params['p60']
                med = params['p50']

                #Plot regression lines of 7, 30 and 100d trends
                if p in [7, 30, 100]:
                    x = pd.to_datetime([params['s'], params['e']])
                    # Convert time difference into numeric units (e.g., days since x_min)
                    t = (x - params['s']).total_seconds() / (24 * 3600)  # convert to days
                    y = params['b'] * t + params['a']
                    ax2.plot(x, y, label=f'{p}d-Trend (PosRate)', linestyle='dashed', color=TREND_COLORS[p])
        
        # Plot annual mean
        ax2.axhline(med, label="Jahresmedian (PosRate)", linewidth=1, color='#404040', linestyle='dotted')

        # Plot RKI incidence if available for spec
        if RKI and spec in df_rki.columns:
            print('Adding RKI data...')
            ax2.bar(df_rki.index, df_rki[spec], label=f'RKI-Inzidenz (Region: {REGION})', width=0.8, color='green')
            # extra text label on the same right axis
            ax2.text(
                1.03, 0.5, f'RKI-Inzidenz {REGION} (pro 100k)', #position at the very right center
                transform=ax2.transAxes,
                rotation=270,
                va='center',
                ha='left',
                color='green'
            )

        # Add to arrow table 
        htmlRef = f'<td style="color:{col}; padding: 5px;">{curRate:.1f}% ({p40:.1f} - {p60:.1f})</td>'
        htmlBody += f'<tr><td style="color:{col}; padding: 5px;">{spec}</td>' + htmlTrends + htmlRef
        
        if spec not in TOPSPECS:
            inclusion = True
            if curRate < MINRATE:
                txtAdd = f' (aktuell exkludiert, da Positivrate < {MINRATE}%)'
                inclusion = False
            elif MINPERCENTILE == 'p60' and curRate <= p60:
                txtAdd = f' (aktuell exkludiert, da Positivrate <= {MINPERCENTILE} von {p60:.1f}%)'
                inclusion = False
            if not inclusion:
                axs[plotCounter].set_title(spec + txtAdd)
                htmlBody += f'<td style="font-size: 9px;">{txtAdd}</td>'
        htmlBody += '</tr>'
        #Keep the handles an labels of the Corona plot (valid for all other)
        if plotCounter == 0:
            # Get handles and labels from both subplots
            handles, labels = axs[0].get_legend_handles_labels()
            handles2, labels2 = ax2.get_legend_handles_labels()
            handles += handles2
            labels += labels2

        plotCounter += 1
    
    
    print(labels)
    # Create a single legend for the entire figure
    topMargin = 0.96
    # specify order of legend labels
    order = [0,1,-1,2,3,4]
    fig.legend([handles[i] for i in order], [labels[i] for i in order], loc='lower center', ncol=len(handles), bbox_to_anchor=(0.5, 1), bbox_transform=fig.transFigure)
    
    plt.subplots_adjust(top=topMargin) #leave more space for legend on top
    plt.tight_layout()
    # --- save to file ---
    outfile = f'{OUTDIR}/KORLAB_RespVir_{PERIOD}_{end}.png'
    plt.savefig(outfile,       # filename (extensions: .png, .pdf, .svg, etc.)
                dpi=80,                      # dots per inch (resolution); default: 80
                bbox_inches='tight')          # trim excess margins
    print(f'Plot saved to {outfile}')
    fp = open(outfile, 'rb')
    img = MIMEImage(fp.read(), name=f'KORLAB_RespVir_{PERIOD}_{end}.png')
    fp.close()
    msg.attach(img)
    print('Generating HTML...')
    htmlBody += '</table>'
    analyteString = '<ul>' # style="font-size: 10px">'
    #for key, uunrs in seenUUNRs.items():
    for inftype, uunrs in uunr_map.items():
        sortedUUNRs = map(str, sorted((map(int, uunrs)))) #sort numerically as ints, then remap to string
        analyteString += f'<li>{inftype}: {", ".join(sortedUUNRs)}</li>'
    analyteString += '</ul>'
    #print(analyteString)

    validString = ', '.join(validResults.keys())
    htmlBody += f'''
        <hr>
        <p>Positivrate auf <font style="color:red;">hohem</font>, <font style="color:orange;">erhöhtem</font>, <font style="color:black;">mittlerem</font>, <font style="color:#aaaaaa;">erniedrigtem</font> oder <font style="color:green;">niedrigem</font> Niveau</p>
        <p>Trend (Änderung in %) über angegebenen Zeitraum: &#8679&#8679: starker Anstieg, &#8679: leichter Anstieg, &#8680: kaum Veränderung, &#8681: leichter Rückgang, &#8681&#8681: starker Rückgang</p>
        <hr>
        <p>Gültige Werte: <i>{validString}</i>. Alle anderen Werte gelten als ungültig und fließen nicht mit in die Berechnungen ein.<br>
        Berücksichtigte Medat-Analyten je Erregertyp:<br>{analyteString}</p>
        <hr>
        <p><i>Diese E-Mail wurde automatisch generiert. Bei Rückfragen senden Sie bitte eine E-Mail an <a href="bioinformatik@laborberlin.com">bioinformatik@laborberlin.com</a></i></p>

        '''  
    print(htmlBody)

    # Send the message via our own SMTP server.
    if SEND_MAIL:
        msg["Subject"] = f'KORLAB: RespVir ({end.strftime("%d.%m.%Y")})'
        msg["From"] = EMAIL_FROM
        msg["To"] = ", ".join(EMAIL_TO)
        msg['Cc'] = ", ".join(EMAIL_CC)
        #msg.add_header("Content-Type", "text/html")
        #msg.set_payload(htmlBody, charset="utf-8")

        
        msg.attach(MIMEText(htmlBody, 'html'))
        s = smtplib.SMTP('localhost')
        s.send_message(msg)
        s.quit()
        print(f'Email sent to {EMAIL_TO}')

if __name__ == "__main__":
    main()
