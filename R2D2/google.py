'''
    Functions for push information to Google Spreadsheet
'''
import glob, os

def init_gspread(json_key, project):
    '''
    This function initialize the utility of google spread 

    Parameters:
        json_key (str): file of json key to access Google API
        projet (str): project name, typically name of upper directory

    Returnes:
        gc (instance): Google API instance
    '''

    import gspread
    from oauth2client.service_account import ServiceAccountCredentials

    scope = ['https://spreadsheets.google.com/feeds',
             'https://www.googleapis.com/auth/drive']

    credentials = ServiceAccountCredentials.from_json_keyfile_name(json_key, scope)
    gc = gspread.authorize(credentials)
    
    return gc
################################################################################
def fetch_URL_gspread(
    json_key=None,
    project=__file__.split('/')[-4],    
    ):
    '''
    fetch corresponding URL for google spreadsheet

    Parameters:
        json_key (str): file of json key to access Google API
        projet (str): project name, typically name of upper directory

    Returnes:
        URL (str): Google spreadsheet URL
    '''
    
    if json_key == None:
        json_key = glob.glob(os.environ['HOME']+'/json/*')[0]
    gid = init_gspread(json_key,project).open(project).id
    
    return 'https://docs.google.com/spreadsheets/d/'+gid
################################################################################

def set_top_line(
    json_key=None,
    project=__file__.split('/')[-4],
    ):
    '''
    This function set top line of google spreadsheet

    Parameters:
        json_key (str): file of json key to access Google API
        projet (str): project name, typically name of upper directory

    Returnes:
        None
    '''
    if json_key == None:
        json_key = glob.glob(os.environ['HOME']+'/json/*')[0]  
    gc = init_gspread(json_key,project)
    wks = gc.open(project).sheet1

    cells = wks.range('A1:T1')
    keys = [ 'Case ID' \
             ,'Mstar' \
             ,'(ix,jx,kx)' \
             ,'xmin [Mm]' \
             ,'xmax [Mm]' \
             ,'ymin' \
             ,'ymax' \
             ,'zmin' \
             ,'zmax' \
             ,'uni' \
             ,'dx [km]' \
             ,'m ray' \
             ,'dtout [s]' \
             ,'dtout_tau [s]' \
             ,'al' \
             ,'RSST' \
             ,'Om' \
             ,'Gemetry' \
             ,'origin'
             ,'update time' \
             ,'Server'
             ]

    for cell, key in zip(cells,keys):
        cell.value = key

    wks.update_cells(cells)
        
################################################################################
def set_cells_gspread(self,
                json_key=None,
                project=__file__.split('/')[-4],
                caseid=None):
    '''
    This function output parameters to 
    Google spread sheet

    Parameters:
        self (R2D2_data): instance of R2D2_data class
        json_key (str): file of json key to access Google API
        projet (str): project name, typically name of upper directory
        caseid (str): caseid
    
    '''
    import datetime
    import numpy as np
    import R2D2

    if json_key == None:
        json_key = glob.glob(os.environ['HOME']+'/json/*')[0]  

    if caseid is None:
        caseid = self.p['datadir'].split('/')[-3]
    
    gc = init_gspread(json_key,project)
    wks = gc.open(project).sheet1
    str_id = str(int(caseid[1:])+1)
    cells = wks.range('A'+str_id+':'+'T'+str_id)

    keys = [caseid]
    keys.append('{:.2f}'.format(self.p['mstar']/R2D2.msun))
    keys.append(str(self.p['ix'])+' '+str(self.p['jx'])+' '+str(self.p['kx']))
    keys.append( '{:6.2f}'.format((self.p['xmin']-self.p['rstar'])*1.e-8))
    keys.append( '{:6.2f}'.format((self.p['xmax']-self.p['rstar'])*1.e-8))
    

    if self.p['geometry'] == 'Cartesian':
        keys.append( '{:6.2f}'.format(self.p['ymin']*1.e-8)+' [Mm]')
        keys.append( '{:6.2f}'.format(self.p['ymax']*1.e-8)+' [Mm]')
        keys.append( '{:6.2f}'.format(self.p['zmin']*1.e-8)+' [Mm]')
        keys.append( '{:6.2f}'.format(self.p['zmax']*1.e-8)+' [Mm]')

    if self.p['geometry'] == 'Spherical':
        pi2rad = 180/np.pi
        keys.append( '{:6.2f}'.format(self.p['ymin']*pi2rad)+' [deg]')
        keys.append( '{:6.2f}'.format(self.p['ymax']*pi2rad)+' [deg]')
        keys.append( '{:6.2f}'.format(self.p['zmin']*pi2rad)+' [deg]')
        keys.append( '{:6.2f}'.format(self.p['zmax']*pi2rad)+' [deg]')

    if self.p['geometry'] == 'YinYang':
        pi2rad = 180/np.pi
        keys.append( '0 [deg]')
        keys.append( '180 [deg]')
        keys.append( '-180 [deg]')
        keys.append( '180 [deg]')
    
    if self.p['ununiform_flag']:
        keys.append('F')
    else:
        keys.append('T')
    
    dx0 = (self.p['x'][1] - self.p['x'][0])*1.e-5
    dx1 = (self.p['x'][self.p['ix']-1] - self.p['x'][self.p['ix']-2])*1.e-5
    keys.append( '{:6.2f}'.format(dx0)+' '+'{:6.2f}'.format(dx1))
    keys.append( self.p['rte'])
    keys.append( '{:6.2f}'.format(self.p['dtout']))
    keys.append( '{:6.2f}'.format(self.p['dtout_tau']))
    keys.append( '{:5.2f}'.format(self.p['potential_alpha']))
    if self.p['xi'].max() == 1.0:
        keys.append('F')
    else:
        keys.append('T')

    keys.append( '{:5.1f}'.format(self.p['omfac']))
    keys.append(self.p['geometry'])
    keys.append(self.p['origin'])
    keys.append(str(datetime.datetime.now()).split('.')[0])
    keys.append(self.p['server'])
    
    for cell, key in zip(cells,keys):
        cell.value = key
    
    wks.update_cells(cells)

