def init(self, datadir, verbose=False, self_old=None):
    '''
    This method reads basic data for calculation setting
    The data is stored in self.p dictionary

    Parameters:
        datadir (str): location of data
        verbose (bool): True shows the information of data
        self_old (R2D2_data): if self_old is not None, datadir is compared with old one and if datadir is same as old one, self is updated with self_old
    '''
    import R2D2
    import numpy as np
    from scipy.io import FortranFile
    import os,sys
    
    # check if the data is already read
    initialize_flag = True
    if self_old is not None:
        if self_old.p['datadir'] == datadir:
            initialize_flag = False
          
    if initialize_flag:  
        self.p = {}
        self.qs = {}
        self.qz = {}
        self.qq = {}
        self.qr = {}
        self.qv = {}
        self.qt = {}
        self.qd = {}
        self.qt_yin = {}
        self.qt_yan = {}
        self.ql = {}
        self.ql_yin = {}
        self.ql_yan = {}
        self.q2 = {}
        self.t = {}
        self.vc = {}
    else:   
        self.__dict__.update(vars(self_old))
            
    self.p['datadir'] = datadir 

    # read basic parameters from param directory
    f = open(self.p['datadir']+"param/nd.dac","r")
    nn = f.read().split()
    nd = int(nn[0]) # current maximum time step
    nd_tau = int(nn[1]) # current maximum time step for tau
    f.close()
    
    self.p['nd'] = nd
    # data/time/tau のファイルの数を数え、nd_tauと比較して大きい方をnd_tauとする
    nd_tau = max(len(os.listdir(datadir+'time/tau')) - 1,nd_tau)
    self.p['nd_tau'] = nd_tau

    ## version check
    R2D2_py_ver = 2.0
    f = open(self.p['datadir']+"param/params.dac","r")
    line = f.readline().split()
    if R2D2_py_ver != float(line[2]):
        print("#######################################################")
        print("#######################################################")
        print("### Current R2D2 Python version is ",R2D2_py_ver,".")
        print("### You use the data from R2D2 version ",float(line[2]),".")
        print("##p# Please use the same version of fortran R2D2.")
        print("#######################################################")
        print("#######################################################")
        sys.exit()
    
    # Basic parameterの読み込み
    line = f.readline()
    while line:
        if line.split()[2] == 'i': # integer
            self.p[line.split()[1]] = int(line.split()[0])
        if line.split()[2] == 'd': # double
            self.p[line.split()[1]] = float(line.split()[0])
        if line.split()[2] == 'c': # character
            self.p[line.split()[1]] = line.split()[0]
        if line.split()[2] == 'l':
            if line.split()[0] == 'F':
                self.p[line.split()[1]] = False
            else:
                self.p[line.split()[1]] = True
                
        line = f.readline()
    f.close()

    # transform endiant for python
    if self.p["swap"] == 0: # little
        self.p["endian"] = "<"
    else:
        self.p["endian"] = ">"

    # MPI information
    f = FortranFile(self.p['datadir']+'param/xyz.dac','r')
    shape = (self.p['ix0']*self.p['jx0']*self.p['kx0'],3)
    self.p['xyz'] = f.read_reals(dtype=np.int32).reshape(shape,order='F')
    ## number of grid in each direction
    self.p["ix"] = self.p["ix0"]*self.p["nx"]
    self.p["jx"] = self.p["jx0"]*self.p["ny"]
    self.p["kx"] = self.p["kx0"]*self.p["nz"]
   
    ixg = self.p["ix"] + 2*self.p["margin"]*(self.p["xdcheck"]-1)
    jxg = self.p["jx"] + 2*self.p["margin"]*(self.p["ydcheck"]-1)
    kxg = self.p["kx"] + 2*self.p["margin"]*(self.p["zdcheck"]-1)
        
    # read background variables
    endian = self.p["endian"]
    dtyp=np.dtype([ \
                    ("head",endian+"i",),\
                    ("x",endian+"d",(ixg,)),\
                    ("y",endian+"d",(jxg,)),\
                    ("z",endian+"d",(kxg,)),\
                    ("pr0",endian+"d",(ixg,)),\
                    ("te0",endian+"d",(ixg,)),\
                    ("ro0",endian+"d",(ixg,)),\
                    ("se0",endian+"d",(ixg,)),\
                    ("en0",endian+"d",(ixg,)),\
                    ("op0",endian+"d",(ixg,)),\
                    ("tu0",endian+"d",(ixg,)),\
                    ("dsedr0",endian+"d",(ixg,)),\
                    ("dtedr0",endian+"d",(ixg,)),\
                    ("dprdro",endian+"d",(ixg,)),\
                    ("dprdse",endian+"d",(ixg,)),\
                    ("dtedro",endian+"d",(ixg,)),\
                    ("dtedse",endian+"d",(ixg,)),\
                    ("dendro",endian+"d",(ixg,)),\
                    ("dendse",endian+"d",(ixg,)),\
                    ("gx",endian+"d",(ixg,)),\
                    ("cp",endian+"d",(ixg,)),\
                    ("fa",endian+"d",(ixg,)),\
                    ("sa",endian+"d",(ixg,)),\
                    ("xi",endian+"d",(ixg,)),\
                    ("tail",endian+"i")\
    ])
    f = open(self.p['datadir']+"param/back.dac",'rb')
    back = np.fromfile(f,dtype=dtyp,count=1)
    f.close()

    # marginも含んだ座標
    self.p['xg'] = back['x'].reshape((ixg),order='F')
    self.p['yg'] = back['y'].reshape((jxg),order='F')
    self.p['zg'] = back['z'].reshape((kxg),order='F')

    # fluxはi+1/2で定義するので、そのための座標を定義
    self.p['x_flux'] = np.zeros(self.p["ix"]+1)
    for i in range(0,self.p["ix"]+1):
        self.p['x_flux'][i] = 0.5*(self.p['xg'][i+self.p['margin']] + self.p['xg'][i+self.p['margin']-1])
    
    self.p['ro0g'] = back['ro0'].reshape((ixg),order='F')
    self.p['se0g'] = back['se0'].reshape((ixg),order='F')
    
    for key in back.dtype.names:
        if back[key].size == ixg:
            self.p[key] = back[key].reshape((ixg),order="F")[self.p["margin"]:ixg-self.p["margin"]]
        elif back[key].size == jxg:
            self.p[key] = back[key].reshape((jxg),order="F")[self.p["margin"]:jxg-self.p["margin"]]
        elif back[key].size == kxg:
            self.p[key] = back[key].reshape((kxg),order="F")[self.p["margin"]:kxg-self.p["margin"]]
            
    self.p["xr"] = self.p["x"]/self.p["rstar"]

    if self.p["geometry"] == 'YinYang':
        
        # original geometry in Yin-Yang grid
        self.p['jx_yy'] = self.p['jx']
        self.p['kx_yy'] = self.p['kx']

        self.p['jxg_yy'] = self.p['jx'] + 2*self.p['margin']
        self.p['kxg_yy'] = self.p['kx'] + 2*self.p['margin']
        
        self.p['y_yy'] = self.p['y'] - 0.5*(self.p['y'][1] - self.p['y'][0])
        self.p['z_yy'] = self.p['z'] - 0.5*(self.p['z'][1] - self.p['z'][0])

        self.p['yg_yy'] = self.p['yg']
        self.p['zg_yy'] = self.p['zg']
        
        # converted spherical geometry
        self.p['jx'] = self.p['jx']*2
        self.p['kx'] = self.p['jx']*2

        dy =   np.pi/self.p['jx']
        dz = 2*np.pi/self.p['kx']

        y = np.zeros(self.p['jx'])
        z = np.zeros(self.p['kx'])
        
        y[0] = 0.5*dy
        z[0] = 0.5*dz - np.pi

        for j in range(1,self.p['jx']):
            y[j] = y[j-1] + dy

        for k in range(1,self.p['kx']):
            z[k] = z[k-1] + dz

        self.p['y'] = y
        self.p['z'] = z
            
    if self.p['zdcheck'] == 2:
        dimension = '3d'
    else:
        dimension = '2d'
        
    ##############################
    # read value information
    if dimension == "3d":
        f = open(self.p['datadir']+"remap/vl/c.dac","r")
        value = f.read().split()
        self.p["m2d_xy"] = int(value[0])
        self.p["m2d_xz"] = int(value[1])
        self.p["m2d_flux"] = int(value[2])
        if self.p['geometry'] == 'YinYang':
            self.p["m2d_spex"] = int(value[3])
        del value[0:3]
        if self.p['geometry'] == 'YinYang':
            del value[0]
        self.p["cl"] = list(map(str.strip,value)) ## strip space from character
        f.close()

        ##############################
        # read mpi information for remap
        npe = self.p["npe"]
        dtyp=np.dtype([ \
                    ("iss",endian+str(npe)+"i4"),\
                    ("iee",endian+str(npe)+"i4"),\
                    ("jss",endian+str(npe)+"i4"),\
                    ("jee",endian+str(npe)+"i4"),\
                    ("iixl",endian+str(npe)+"i4"),\
                    ("jjxl",endian+str(npe)+"i4"),\
                    ("np_ijr",endian+str(self.p["ixr"]*self.p["jxr"])+"i4"),\
                    ("ir",endian+str(npe)+"i4"),\
                    ("jr",endian+str(npe)+"i4"),\
                    ("i2ir",endian+str(ixg)+"i4"),\
                    ("j2jr",endian+str(jxg)+"i4"),\
        ])
        
        f = open(datadir+"remap/remap_info.dac",'rb')
        mpi = np.fromfile(f,dtype=dtyp,count=1)    
        f.close()
        
        for key in mpi.dtype.names:
            if key == "np_ijr":
                self.p[key] = mpi[key].reshape((self.p["ixr"],self.p["jxr"]),order="F")
            else:
                self.p[key] = mpi[key].reshape((mpi[key].size),order="F")
    
        self.p["i2ir"] = self.p["i2ir"][self.p["margin"]:ixg-self.p["margin"]]
        self.p["j2jr"] = self.p["j2jr"][self.p["margin"]:jxg-self.p["margin"]]
        
        self.p["iss"] = self.p["iss"] - 1
        self.p["iee"] = self.p["iee"] - 1
        self.p["jss"] = self.p["jss"] - 1
        self.p["jee"] = self.p["jee"] - 1

        if os.path.isdir(self.p['datadir']+'slice'):
            f = open(self.p['datadir']+"slice/params.dac","r")
            line = f.readline()
            while line:
                self.p[line.split()[1]] = int(line.split()[0])
                line = f.readline()

            f.close()

            dtype = np.dtype([ \
                ('x_slice',self.p['endian']+str(self.p['nx_slice'])+'d'),\
                ('y_slice',self.p['endian']+str(self.p['ny_slice'])+'d'),\
                ('z_slice',self.p['endian']+str(self.p['nz_slice'])+'d'),\
            ])
            f = open(self.p['datadir']+'slice/slice.dac','rb')
            slice = np.fromfile(f,dtype=dtype)

            self.p['x_slice'] = slice['x_slice'].reshape((self.p['nx_slice']),order='F')
            self.p['y_slice'] = slice['y_slice'].reshape((self.p['ny_slice']),order='F')
            self.p['z_slice'] = slice['z_slice'].reshape(self.p['nz_slice'],order='F')
            
            f.close()
    # read order data
    srcdir = datadir[:-5]+'src/all/'
    if os.path.exists(srcdir+'info.txt'):
        f = open(srcdir+'info.txt')
        orders = f.read().split(',')
        self.p['order_1D'] = int(orders[0])
        self.p['order_2D'] = int(orders[1])
        self.p['order_3D'] = int(orders[2])
    else:
        order_3D = 1
        
    # read equation of state
    eosdir = datadir[:-5]+'input_data/'
    if os.path.exists(eosdir+'eos_table_sero.npz'):
        eos_d = np.load(eosdir+'eos_table_sero.npz')
        self.p['log_ro_e'] = eos_d['ro'] # density is defined in logarithmic scale
        self.p['se_e'] = eos_d['se']
        self.p['ix_e'] = len(self.p['log_ro_e'])
        self.p['jx_e'] = len(self.p['se_e'])
        
        self.p['log_pr_e'] = np.log(eos_d['pr']+ 1.e-200)
        self.p['log_en_e'] = np.log(eos_d['en']+ 1.e-200)
        self.p['log_te_e'] = np.log(eos_d['te']+ 1.e-200)
        self.p['log_op_e'] = np.log(eos_d['op']+ 1.e-200)
        
        self.p['dlogro_e'] = self.p['log_ro_e'][1] - self.p['log_ro_e'][0]
        self.p['dse_e'] = self.p['se_e'][1] - self.p['se_e'][0]
                
    # read original data
    if os.path.exists(datadir+'cont_log.txt'):
        f = open(datadir+'cont_log.txt')
        self.p['origin'] = f.readlines()[6][-11:-7]
    else:
        self.p['origin'] = 'N/A'
            
    if verbose:
        R2D2.util.show_information(self)
        
##############################
def read_qq_select(self,xs,n,silent=False):
    '''
    This method reads 2D slice data at a selected height.
    The data is stored in self.qs dictionary

    Parameters:
        xs (float): a selected height for data
        n (int): a setected time step for data
        silent (bool): True suppresses a message of store
        
    '''

    import numpy as np
    i0 = np.argmin(np.abs(self.p["x"]-xs))
    ir0 = self.p["i2ir"][i0]
    mtype = self.p["mtype"]
    iixl = self.p["iixl"]
    jjxl = self.p["jjxl"]
    jx = self.p["jx"]
    kx = self.p["kx"]
    iss = self.p["iss"]
    jss = self.p["jss"]
    jee = self.p["jee"]
    order_3D = self.p["order_3D"]
    
    ### Only when memory is not allocated 
    ### and the size of array is different
    ### memory is allocated
    memflag = True
    if 'ro' in self.qs:
        memflag = not self.qs['ro'].shape == (jx,kx)
    if 'ro' not in self.qs or memflag:
        print('memory is newly allocated')
        self.qs["ro"] = np.zeros((jx,kx))
        self.qs["vx"] = np.zeros((jx,kx))
        self.qs["vy"] = np.zeros((jx,kx))
        self.qs["vz"] = np.zeros((jx,kx))
        self.qs["bx"] = np.zeros((jx,kx))
        self.qs["by"] = np.zeros((jx,kx))
        self.qs["bz"] = np.zeros((jx,kx))
        self.qs["se"] = np.zeros((jx,kx))
        self.qs["ph"] = np.zeros((jx,kx))
        self.qs["pr"] = np.zeros((jx,kx))
        self.qs["te"] = np.zeros((jx,kx))
        self.qs["op"] = np.zeros((jx,kx))
    for jr0 in range(1,self.p["jxr"]+1):
        np0 = self.p["np_ijr"][ir0-1,jr0-1]
 
        if jr0 == self.p['jr'][np0]:
            dtype=np.dtype([ \
                    ("qq",self.p["endian"]+str(mtype*iixl[np0]*jjxl[np0]*kx)+"f"),\
                    ("pr",self.p["endian"]+str(iixl[np0]*jjxl[np0]*kx)+"f"),\
                    ("te",self.p["endian"]+str(iixl[np0]*jjxl[np0]*kx)+"f"),\
                    ("op",self.p["endian"]+str(iixl[np0]*jjxl[np0]*kx)+"f"),\
            ])        
            cnou = '{0:05d}'.format(np0//1000)
            cno  = '{0:08d}'.format(np0)
            f = open(self.p['datadir']+"remap/qq/"+cnou+"/"+cno+"/qq.dac."+'{0:08d}'.format(n)+"."+'{0:08d}'.format(np0),'rb')
            qqq = np.fromfile(f,dtype=dtype,count=1)
            index = [iixl[np0], jjxl[np0], kx, mtype]
            self.qs["ro"][jss[np0]:jee[np0]+1,:] = qqq["qq"].reshape((iixl[np0],jjxl[np0],kx,mtype),order="F")[i0-iss[np0],:,:,0] 
            self.qs["vx"][jss[np0]:jee[np0]+1,:] = qqq["qq"].reshape((iixl[np0],jjxl[np0],kx,mtype),order="F")[i0-iss[np0],:,:,1] 
            self.qs["vy"][jss[np0]:jee[np0]+1,:] = qqq["qq"].reshape((iixl[np0],jjxl[np0],kx,mtype),order="F")[i0-iss[np0],:,:,2] 
            self.qs["vz"][jss[np0]:jee[np0]+1,:] = qqq["qq"].reshape((iixl[np0],jjxl[np0],kx,mtype),order="F")[i0-iss[np0],:,:,3] 
            self.qs["bx"][jss[np0]:jee[np0]+1,:] = qqq["qq"].reshape((iixl[np0],jjxl[np0],kx,mtype),order="F")[i0-iss[np0],:,:,4] 
            self.qs["by"][jss[np0]:jee[np0]+1,:] = qqq["qq"].reshape((iixl[np0],jjxl[np0],kx,mtype),order="F")[i0-iss[np0],:,:,5] 
            self.qs["bz"][jss[np0]:jee[np0]+1,:] = qqq["qq"].reshape((iixl[np0],jjxl[np0],kx,mtype),order="F")[i0-iss[np0],:,:,6] 
            self.qs["se"][jss[np0]:jee[np0]+1,:] = qqq["qq"].reshape((iixl[np0],jjxl[np0],kx,mtype),order="F")[i0-iss[np0],:,:,7] 
            self.qs["ph"][jss[np0]:jee[np0]+1,:] = qqq["qq"].reshape((iixl[np0],jjxl[np0],kx,mtype),order="F")[i0-iss[np0],:,:,8]
   
            self.qs["pr"][jss[np0]:jee[np0]+1,:] = qqq["pr"].reshape((iixl[np0],jjxl[np0],kx),order="F")[i0-iss[np0],:,:]
            self.qs["te"][jss[np0]:jee[np0]+1,:] = qqq["te"].reshape((iixl[np0],jjxl[np0],kx),order="F")[i0-iss[np0],:,:]
            self.qs["op"][jss[np0]:jee[np0]+1,:] = qqq["op"].reshape((iixl[np0],jjxl[np0],kx),order="F")[i0-iss[np0],:,:]
            info = {}
            info['xs'] = self.p['x'][i0]
            self.qs['info'] = info
            f.close()        

    if not silent :
        print('### variables are stored in self.qs ###')

##############################
def read_qq_select_z(self,zs,n,silent=False):
    '''
    This method reads 2D slice data at a selected height.
    The data is stored in self.qs dictionary

    Parameters:
        zs (float): a selected z for data
        n (int): a selected time step for data
        silent (bool): True suppresses a message of store
        
    '''

    import numpy as np
    k0 = np.argmin(np.abs(self.p["z"]-zs))
    mtype = self.p["mtype"]
    iixl = self.p["iixl"]
    jjxl = self.p["jjxl"]
    ix = self.p["ix"]
    jx = self.p["jx"]
    kx = self.p["kx"]
    iss = self.p["iss"]
    iee = self.p["iee"]
    jss = self.p["jss"]
    jee = self.p["jee"]
    order_3D = self.p["order_3D"]
    
    ### Only when memory is not allocated 
    ### and the size of array is different
    ### memory is allocated
    memflag = True
    if 'ro' in self.qz:
        memflag = not self.qz['ro'].shape == (ix,jx)
    if 'ro' not in self.qz or memflag:
        print('memory is newly allocated')
        self.qz["ro"] = np.zeros((ix,jx))
        self.qz["vx"] = np.zeros((ix,jx))
        self.qz["vy"] = np.zeros((ix,jx))
        self.qz["vz"] = np.zeros((ix,jx))
        self.qz["bx"] = np.zeros((ix,jx))
        self.qz["by"] = np.zeros((ix,jx))
        self.qz["bz"] = np.zeros((ix,jx))
        self.qz["se"] = np.zeros((ix,jx))
        self.qz["ph"] = np.zeros((ix,jx))
        self.qz["pr"] = np.zeros((ix,jx))
        self.qz["te"] = np.zeros((ix,jx))
        self.qz["op"] = np.zeros((ix,jx))

    for ir0 in range(1,self.p["ixr"]+1):
        print(ir0,self.p['ixr'])
        for jr0 in range(1,self.p["jxr"]+1):
            np0 = self.p["np_ijr"][ir0-1,jr0-1]
            dtyp=np.dtype([ \
                            ("qq",self.p["endian"]+str(mtype*iixl[np0]*jjxl[np0]*kx)+"f"),\
                            ("pr",self.p["endian"]+str(iixl[np0]*jjxl[np0]*kx)+"f"),\
                            ("te",self.p["endian"]+str(iixl[np0]*jjxl[np0]*kx)+"f"),\
                            ("op",self.p["endian"]+str(iixl[np0]*jjxl[np0]*kx)+"f"),\
            ])
            f = open(self.p['datadir']+"remap/qq/qq.dac."+'{0:08d}'.format(n)+"."+'{0:08d}'.format(np0),'rb')
            qqq = np.fromfile(f,dtype=dtyp,count=1)
            if(order_3D == 4):
                index = [iixl[np0], jjxl[np0], kx, mtype]
                self.qz["ro"][iss[np0]:iee[np0]+1,jss[np0]:jee[np0]+1] = qqq["qq"].reshape((iixl[np0], jjxl[np0], kx, mtype), order="F")[:,:,k0,0] 
                self.qz["vx"][iss[np0]:iee[np0]+1,jss[np0]:jee[np0]+1] = qqq["qq"].reshape((iixl[np0], jjxl[np0], kx, mtype), order="F")[:,:,k0,1] 
                self.qz["vy"][iss[np0]:iee[np0]+1,jss[np0]:jee[np0]+1] = qqq["qq"].reshape((iixl[np0], jjxl[np0], kx, mtype), order="F")[:,:,k0,2]
                self.qz["vz"][iss[np0]:iee[np0]+1,jss[np0]:jee[np0]+1] = qqq["qq"].reshape((iixl[np0], jjxl[np0], kx, mtype), order="F")[:,:,k0,3]
                self.qz["bx"][iss[np0]:iee[np0]+1,jss[np0]:jee[np0]+1] = qqq["qq"].reshape((iixl[np0], jjxl[np0], kx, mtype), order="F")[:,:,k0,4]
                self.qz["by"][iss[np0]:iee[np0]+1,jss[np0]:jee[np0]+1] = qqq["qq"].reshape((iixl[np0], jjxl[np0], kx, mtype), order="F")[:,:,k0,5]
                self.qz["bz"][iss[np0]:iee[np0]+1,jss[np0]:jee[np0]+1] = qqq["qq"].reshape((iixl[np0], jjxl[np0], kx, mtype), order="F")[:,:,k0,6]
                self.qz["se"][iss[np0]:iee[np0]+1,jss[np0]:jee[np0]+1] = qqq["qq"].reshape((iixl[np0], jjxl[np0], kx, mtype), order="F")[:,:,k0,7]
                self.qz["ph"][iss[np0]:iee[np0]+1,jss[np0]:jee[np0]+1] = qqq["qq"].reshape((iixl[np0], jjxl[np0], kx, mtype), order="F")[:,:,k0,8]
            
                self.qz["pr"][iss[np0]:iee[np0]+1,jss[np0]:jee[np0]+1] = qqq["pr"].reshape((iixl[np0],jjxl[np0],kx),order="F")[:,:,k0]
                self.qz["te"][iss[np0]:iee[np0]+1,jss[np0]:jee[np0]+1] = qqq["te"].reshape((iixl[np0],jjxl[np0],kx),order="F")[:,:,k0]
                self.qz["op"][iss[np0]:iee[np0]+1,jss[np0]:jee[np0]+1] = qqq["op"].reshape((iixl[np0],jjxl[np0],kx),order="F")[:,:,k0]

            f.close()
                
        info = {}
        info['zs'] = self.p['z'][k0]
        self.qz['info'] = info
        
    if not silent :
        print('### variables are stored in self.qz ###')
        
##############################
def read_qq(self,n,value,silent=False):
    '''
    This method reads 3D full data
    The data is stored in self.qq dictionary

    Parameters:
        n (int): a selected time step for data
        value (char): kind of value
            "ro": density
            "vx", "vy", "vz": velocity
            "bx", "by", "bz": magnetic field
            "se": entropy
            "ph": div B cleaning
            "te": temperature
            "op": Opacity
        silent (bool): True suppresses a message of store
    '''
    import numpy as np
    
    mtype = self.p["mtype"]
    iixl = self.p["iixl"]
    jjxl = self.p["jjxl"]
    kx = self.p["kx"]
    iss = self.p["iss"]
    iee = self.p["iee"]
    jss = self.p["jss"]
    jee = self.p["jee"]
    ix = self.p["ix"]
    jx = self.p["jx"]
    kx = self.p["kx"]
    order_3D = self.p["order_3D"]

    values = ['ro','vx','vy','vz','bx','by','bz','se','ph','pr','te','op']

    if type(value) == str:
        values_input = [value]
    if type(value) == list:
        values_input = value

    for value in values_input:
        if value not in values:
            print('######')
            print('value =',value)
            print('value should be one of ',values)
            print('return')
            return
    
    ### Only when memory is not allocated 
    ### and the size of array is different
    ### memory is allocated
    for value in values_input:
        memflag = True    
        if value in self.qq:
            memflag = not self.qq[value].shape == (ix,jx,kx)
        if  value not in self.qq or memflag:
            print('memory is newly allocated')
            self.qq[value] = np.zeros((ix,jx,kx),dtype=np.float32)

    for ir0 in range(1,self.p["ixr"]+1):
        for jr0 in range(1,self.p["jxr"]+1):
            np0 = self.p["np_ijr"][ir0-1,jr0-1]
            if ir0 == self.p['ir'][np0] and jr0 == self.p['jr'][np0]:
                dtype=np.dtype([ \
                    ("qq",self.p["endian"]+str(mtype*iixl[np0]*jjxl[np0]*kx)+"f"),\
                    ("pr",self.p["endian"]+str(iixl[np0]*jjxl[np0]*kx)+"f"),\
                    ("te",self.p["endian"]+str(iixl[np0]*jjxl[np0]*kx)+"f"),\
                    ("op",self.p["endian"]+str(iixl[np0]*jjxl[np0]*kx)+"f"),\
                ])

                cnou = '{0:05d}'.format(np0//1000)
                cno  = '{0:08d}'.format(np0)
                f = open(self.p['datadir']+"remap/qq/"+cnou+"/"+cno+"/qq.dac."+'{0:08d}'.format(n)+"."+'{0:08d}'.format(np0),'rb')            
                #f = open(self.p['datadir']+"remap/qq/qq.dac."+'{0:08d}'.format(n)+"."+'{0:08d}'.format(np0),'rb')
                qqq = np.fromfile(f,dtype=dtype,count=1)

                for value in values_input:
                    if value in values[:-3]:
                        m = values.index(value)
                        self.qq[value][iss[np0]:iee[np0]+1,jss[np0]:jee[np0]+1,:] = qqq["qq"].reshape((iixl[np0],jjxl[np0],kx,mtype),order="F")[:,:,:,m]
                    else:
                        self.qq[value][iss[np0]:iee[np0]+1,jss[np0]:jee[np0]+1,:] = qqq[value].reshape((iixl[np0],jjxl[np0],kx),order="F")[:,:,:]
                f.close()

    if not silent :
        print('### variables are stored in self.qq ###')
        
##############################
def read_qq_restricted(self,n,value,x0,x1,y0,y1,z0,z1,silent=False):
    '''
    This method reads 3D restricted-area data
    The data is stored in self.qr dictionary

    Parameters:
        n (int): a selected time step for data
        value (char): kind of value
        "ro": density
        "vx", "vy", "vz": velocity
        "bx", "by", "bz": magnetic field
        "se": entropy
        "ph": div B cleaning
        "te": temperature
        "op": Opacity        
        x0, y0, z0 (float): minimum x, y, z
        x1, y1, z1 (float): maximum x, y, z
        silent (bool): True suppresses a message of store
    '''
    import numpy as np
    
    mtype = self.p["mtype"]
    iixl, jjxl = self.p["iixl"], self.p["jjxl"]
    iss, iee = self.p["iss"], self.p["iee"]
    jss, jee = self.p["jss"], self.p["jee"]
    ix,jx,kx = self.p["ix"],self.p["jx"],self.p["kx"]
    x,y,z = self.p['x'],self.p['y'],self.p['z']
    
    i0, i1 = np.argmin(abs(x-x0)), np.argmin(abs(x-x1))
    j0, j1 = np.argmin(abs(y-y0)), np.argmin(abs(y-y1))
    k0, k1 = np.argmin(abs(z-z0)), np.argmin(abs(z-z1))
    ixr = i1 - i0 + 1
    jxr = j1 - j0 + 1
    kxr = k1 - k0 + 1
    
    values = ['ro','vx','vy','vz','bx','by','bz','se','ph','pr','te','op']

    if type(value) == str:
        values_input = [value]
    if type(value) == list:
        values_input = value
    
    for value in values_input:
        self.qr[value] = np.zeros((ixr,jxr,kxr),dtype=np.float32)
        
    for ir0 in range(1,self.p["ixr"]+1):        
        for jr0 in range(1,self.p["jxr"]+1):
            np0 = self.p["np_ijr"][ir0-1,jr0-1]
            
            if(not (iss[np0] > i1 or iee[np0] < i0 or jss[np0] > j1 or jee[np0] < j0) ):
                
                dtyp=np.dtype([ \
                    ("qq",self.p["endian"]+str(mtype*iixl[np0]*jjxl[np0]*kx)+"f"),\
                    ("pr",self.p["endian"]+str(iixl[np0]*jjxl[np0]*kx)+"f"),\
                    ("te",self.p["endian"]+str(iixl[np0]*jjxl[np0]*kx)+"f"),\
                    ("op",self.p["endian"]+str(iixl[np0]*jjxl[np0]*kx)+"f"),\
                ])

                cnou = '{0:05d}'.format(np0//1000)
                cno  = '{0:08d}'.format(np0)
                f = open(self.p['datadir']+"remap/qq/"+cnou+"/"+cno+"/qq.dac."+'{0:08d}'.format(n)+"."+'{0:08d}'.format(np0),'rb')      
                #f = open(self.p['datadir']+"remap/qq/qq.dac."+'{0:08d}'.format(n)+"."+'{0:08d}'.format(np0),'rb')
                qqq = np.fromfile(f,dtype=dtyp,count=1)

                for value in values_input:
                    if value in values[:-3]:
                        m = values.index(value)
                        isrt_rcv = max([0  ,iss[np0]-i0  ])
                        iend_rcv = min([ixr,iee[np0]-i0+1])
                        jsrt_rcv = max([0  ,jss[np0]-j0  ])
                        jend_rcv = min([jxr,jee[np0]-j0+1])
                        
                        isrt_snd = isrt_rcv - (iss[np0]-i0)
                        iend_snd = isrt_snd + (iend_rcv - isrt_rcv)
                        jsrt_snd = jsrt_rcv - (jss[np0]-j0)
                        jend_snd = jsrt_snd + (jend_rcv - jsrt_rcv)

                        self.qr[value][isrt_rcv:iend_rcv,jsrt_rcv:jend_rcv,:] = \
                            qqq["qq"].reshape((iixl[np0],jjxl[np0],kx,mtype),order="F")[isrt_snd:iend_snd,jsrt_snd:jend_snd,k0:k1+1,m]
                    else:
                        self.qr[value][isrt_rcv:iend_rcv,jsrt_rcv:jend_rcv,:] = \
                            qqq[value].reshape((iixl[np0],jjxl[np0],kx),order="F")[isrt_snd:iend_snd,jsrt_snd:jend_snd,k0:k1+1]                         
                f.close()

    if not silent :
        print('### variales are stored in self.qr ###')

##############################
def read_qq_all(self,n,silent=False):
    '''
    This method reads 3D full data all
    The data is stored in self.qq dictionary

    Parameters:
        n (int): a selected time step for data
        silent (bool): True suppresses a message of store
    '''
    value = ['ro','vx','vy','vz','bx','by','bz','se','ph','pr','te','op']
    read_qq(self,n,value,silent=silent)


##############################
def read_qq_tau(self,n,silent=False):
    '''
    This method reads 2D data at certain optical depths.
    The data is stored in self.qt dictionary.
    In this version the selected optical depth is 1, 0.1, and 0.01

    Parameters:
        n (int): a selected time step for data
        silent (bool): True suppresses a message of store
    '''
    import numpy as np

    if self.p['geometry'] == 'YinYang':
        files = ['_yin','_yan']
        qts = [self.qt_yin,self.qt_yan]
    else:
        files=['']
        qts = [self.qt]

    m_tu = self.p['m_tu']
    m_in = self.p['m_in']

    if(self.p['geometry'] == 'YinYang'):
        jx = self.p['jx_yy']
        kx = self.p['kx_yy']
    else:
        jx = self.p['jx']
        kx = self.p['kx']
        
    for file,qt in zip(files,qts):
        f = open(self.p['datadir']+"tau/qq"+file+".dac."+'{0:08d}'.format(n),"rb")
        qq_in0 = np.fromfile(f,self.p["endian"]+'f',self.p["m_tu"]*self.p["m_in"]*jx*kx)
        f.close()
        
        qt["in"] = np.reshape(qq_in0,(m_tu,m_in,jx,kx),order="F")[0,0,:,:]
        qt["ro"] = np.reshape(qq_in0,(m_tu,m_in,jx,kx),order="F")[0,1,:,:]
        qt["se"] = np.reshape(qq_in0,(m_tu,m_in,jx,kx),order="F")[0,2,:,:]
        qt["pr"] = np.reshape(qq_in0,(m_tu,m_in,jx,kx),order="F")[0,3,:,:]
        qt["te"] = np.reshape(qq_in0,(m_tu,m_in,jx,kx),order="F")[0,4,:,:]
        qt["vx"] = np.reshape(qq_in0,(m_tu,m_in,jx,kx),order="F")[0,5,:,:]
        qt["vy"] = np.reshape(qq_in0,(m_tu,m_in,jx,kx),order="F")[0,6,:,:]
        qt["vz"] = np.reshape(qq_in0,(m_tu,m_in,jx,kx),order="F")[0,7,:,:]
        qt["bx"] = np.reshape(qq_in0,(m_tu,m_in,jx,kx),order="F")[0,8,:,:]
        qt["by"] = np.reshape(qq_in0,(m_tu,m_in,jx,kx),order="F")[0,9,:,:]
        qt["bz"] = np.reshape(qq_in0,(m_tu,m_in,jx,kx),order="F")[0,10,:,:]
        qt["he"] = np.reshape(qq_in0,(m_tu,m_in,jx,kx),order="F")[0,11,:,:]
        qt["fr"] = np.reshape(qq_in0,(m_tu,m_in,jx,kx),order="F")[0,12,:,:]

        qt["in01"] = np.reshape(qq_in0,(m_tu,m_in,jx,kx),order="F")[1,0,:,:]
        qt["ro01"] = np.reshape(qq_in0,(m_tu,m_in,jx,kx),order="F")[1,1,:,:]
        qt["se01"] = np.reshape(qq_in0,(m_tu,m_in,jx,kx),order="F")[1,2,:,:]
        qt["pr01"] = np.reshape(qq_in0,(m_tu,m_in,jx,kx),order="F")[1,3,:,:]
        qt["te01"] = np.reshape(qq_in0,(m_tu,m_in,jx,kx),order="F")[1,4,:,:]
        qt["vx01"] = np.reshape(qq_in0,(m_tu,m_in,jx,kx),order="F")[1,5,:,:]
        qt["vy01"] = np.reshape(qq_in0,(m_tu,m_in,jx,kx),order="F")[1,6,:,:]
        qt["vz01"] = np.reshape(qq_in0,(m_tu,m_in,jx,kx),order="F")[1,7,:,:]
        qt["bx01"] = np.reshape(qq_in0,(m_tu,m_in,jx,kx),order="F")[1,8,:,:]
        qt["by01"] = np.reshape(qq_in0,(m_tu,m_in,jx,kx),order="F")[1,9,:,:]
        qt["bz01"] = np.reshape(qq_in0,(m_tu,m_in,jx,kx),order="F")[1,10,:,:]
        qt["he01"] = np.reshape(qq_in0,(m_tu,m_in,jx,kx),order="F")[1,11,:,:]
        qt["fr01"] = np.reshape(qq_in0,(m_tu,m_in,jx,kx),order="F")[1,12,:,:]

        qt["in001"] = np.reshape(qq_in0,(m_tu,m_in,jx,kx),order="F")[2,0,:,:]
        qt["ro001"] = np.reshape(qq_in0,(m_tu,m_in,jx,kx),order="F")[2,1,:,:]
        qt["se001"] = np.reshape(qq_in0,(m_tu,m_in,jx,kx),order="F")[2,2,:,:]
        qt["pr001"] = np.reshape(qq_in0,(m_tu,m_in,jx,kx),order="F")[2,3,:,:]
        qt["te001"] = np.reshape(qq_in0,(m_tu,m_in,jx,kx),order="F")[2,4,:,:]
        qt["vx001"] = np.reshape(qq_in0,(m_tu,m_in,jx,kx),order="F")[2,5,:,:]
        qt["vy001"] = np.reshape(qq_in0,(m_tu,m_in,jx,kx),order="F")[2,6,:,:]
        qt["vz001"] = np.reshape(qq_in0,(m_tu,m_in,jx,kx),order="F")[2,7,:,:]
        qt["bx001"] = np.reshape(qq_in0,(m_tu,m_in,jx,kx),order="F")[2,8,:,:]
        qt["by001"] = np.reshape(qq_in0,(m_tu,m_in,jx,kx),order="F")[2,9,:,:]
        qt["bz001"] = np.reshape(qq_in0,(m_tu,m_in,jx,kx),order="F")[2,10,:,:]
        qt["he001"] = np.reshape(qq_in0,(m_tu,m_in,jx,kx),order="F")[2,11,:,:]
        qt["fr001"] = np.reshape(qq_in0,(m_tu,m_in,jx,kx),order="F")[2,12,:,:]

    if not silent :
        if self.p['geometry'] == 'YinYang':
            print('### variables are stored in self.qt_yin and self.qt_yan ###')
        else:
            print('### variables are stored in self.qt ###')
                
##############################
def read_time(self,n,tau=False,silent=True):
    '''
    This method reads time at a selected time step
    The data is stored in self.t

    Parameters:
        n (int): a selected time step for data
        tau (bool): if True time for optical depth
        silent (bool): True suppresses a message of store

    '''

    import numpy as np

    if tau:
        f = open(self.p['datadir']+"time/tau/t.dac."+'{0:08d}'.format(n),"rb")
        self.t = np.fromfile(f,self.p['endian']+'d',1).reshape((1),order='F')[0]            
        f.close()    
    else:
        f = open(self.p['datadir']+"time/mhd/t.dac."+'{0:08d}'.format(n),"rb")
        self.t = np.fromfile(f,self.p['endian']+'d',1).reshape((1),order='F')[0]            
        f.close()    
    
    if not silent :
        print('### variales are stored in self.t ###')

    return self.t
##############################
def read_vc(self,n,silent=False):
    '''
    This method reads on the fly analysis data from fortran.
    The data is stored in self.vc dictionary

    Parameters:
        n (int): a selected time step for data
        silent (bool): True suppresses a message of store
    '''

    import numpy as np

    # read xy plane data
    f = open(self.p['datadir']+"remap/vl/vl_xy.dac."+'{0:08d}'.format(n),"rb")
    vl = np.fromfile(f,self.p["endian"]+'f',self.p['m2d_xy']*self.p['ix']*self.p['jx']) \
           .reshape((self.p['ix'],self.p['jx'],self.p['m2d_xy']),order="F")
    f.close()

    for m in range(self.p["m2d_xy"]):
        self.vc[self.p["cl"][m]] = vl[:,:,m]

    # read xz plane data
    f = open(self.p['datadir']+"remap/vl/vl_xz.dac."+'{0:08d}'.format(n),"rb")
    vl = np.fromfile(f,self.p["endian"]+'f',self.p['m2d_xz']*self.p['ix']*self.p['kx']) \
           .reshape((self.p['ix'],self.p['kx'],self.p['m2d_xz']),order="F")
    f.close()

    for m in range(self.p["m2d_xz"]):
        self.vc[self.p["cl"][m + self.p['m2d_xy']]] = vl[:,:,m]
        
    # read flux related value
    f = open(self.p['datadir']+"remap/vl/vl_flux.dac."+'{0:08d}'.format(n),"rb")
    vl = np.fromfile(f,self.p["endian"]+'f',self.p['m2d_flux']*(self.p['ix']+1)*self.p['jx']) \
           .reshape((self.p['ix']+1,self.p['jx'],self.p['m2d_flux']),order="F")
    f.close()

    for m in range(self.p["m2d_flux"]):
        self.vc[self.p["cl"][m+self.p["m2d_xy"]+self.p['m2d_xz']]] = vl[:,:,m]
    
    # read spectra
    if self.p['geometry'] == 'YinYang':
        f = open(self.p['datadir']+"remap/vl/vl_spex.dac."+'{0:08d}'.format(n),"rb")
        vl = np.fromfile(f,self.p["endian"]+'f',self.p['m2d_spex']*self.p['ix']*self.p['kx']//4) \
               .reshape((self.p['ix'],self.p['kx']//4,self.p['m2d_spex']),order="F")
        f.close()

        for m in range(self.p["m2d_spex"]):
            self.vc[self.p["cl"][m+self.p["m2d_xy"]+self.p['m2d_xz']+self.p['m2d_flux']]] = vl[:,:,m]
            
    if not silent :
        print('### variables are stored in self.vc ###')

##############################
def read_qq_check(self,n,silent=False,end_step=False):
    '''
    This method reads 3D full data for checkpoint
    The data is stored in self.qc dictionary

    Parameters:
        n (int): a selected time step for data
        silent (bool): True suppresses a message of store
        end_step (bool): If true, checkpoint of end step is read.
    '''

    import numpy as np

    mtype = self.p['mtype']
    ix = self.p['ix']
    jx = self.p['jx']
    kx = self.p['kx']
    margin = self.p['margin']

    ixg = ix + 2*margin*(self.p['xdcheck']-1)
    jxg = jx + 2*margin*(self.p['ydcheck']-1)
    kxg = kx + 2*margin*(self.p['zdcheck']-1)

    order_3D = self.p["order_3D"]
    
    step = '{0:08d}'.format(n)
    if end_step:
        if np.mod(self.p['nd'],2) == 0:
            step = 'e'
        if np.mod(self.p['nd'],2) == 1:
            step = 'o'

    f = open(self.p['datadir']+"qq/qq.dac."+step,'rb')
    if(order_3D == 1):
        self.qc = np.fromfile(f,self.p['endian']+'d',mtype*ixg*jxg*kxg).reshape((mtype,ixg,jxg,kxg),order="F")    
    if(order_3D == 2):
        self.qc = np.fromfile(f,self.p['endian']+'d',mtype*ixg*jxg*kxg).reshape((ixg,mtype,jxg,kxg),order="F")    
    if(order_3D == 3):
        self.qc = np.fromfile(f,self.p['endian']+'d',mtype*ixg*jxg*kxg).reshape((ixg,jxg,mtype,kxg),order="F")    
    if(order_3D == 4):
        self.qc = np.fromfile(f,self.p['endian']+'d',mtype*ixg*jxg*kxg).reshape((ixg,jxg,kxg,mtype),order="F")    

    f.close()
    
    if not silent :
        print('### variables are stored in self.qc ###')
        
##############################
def read_qq_rt(self,n,silent=False):
    '''
    This method reads 3D radiation quantities
    The data is stored in self.qd dictionary

    Parameters:
        n (int): a selected time step for data
        silent (bool): True suppresses a message of store
    '''

    import numpy as np

    mtype = self.p['mtype']
    ix = self.p['ix']
    jx = self.p['jx']
    kx = self.p['kx']
    npe = self.p['npe']
    margin = self.p['margin']
    endian = self.p['endian']
    mrad = 3
    nx = self.p['nx']
    ny = self.p['ny']
    nz = self.p['nz']
    nxg = nx + 2*self.p["margin"]*(self.p["xdcheck"]-1)
    nyg = ny + 2*self.p["margin"]*(self.p["ydcheck"]-1)
    nzg = nz + 2*self.p["margin"]*(self.p["zdcheck"]-1)    

    order_3D = self.p["order_3D"]
    
    step = '{0:08d}'.format(n)

    self.qd['in']        = np.zeros((mrad,2,2,2,ix,jx,kx),dtype=np.float32)
    self.qd['qrad']      = np.zeros((ix,jx,kx),dtype=np.float32)
    self.qd['qrad_1ray'] = np.zeros((ix,jx,kx),dtype=np.float32)

    for np0 in range(npe):
        cnou = '{0:05d}'.format(np0//1000)
        cno  = '{0:08d}'.format(np0)
        print(np0)
        f = open(self.p['datadir']+"qq/"+cnou+"/"+cno+"/rt.dac."+'{0:08d}'.format(n)+"."+'{0:08d}'.format(np0),'rb')
        dtype = np.dtype([ \
            ('in'       ,endian+str(mrad*2*2*2*nxg*nxg*nxg)+'d'),\
            ('qrad'     ,endian+str(           nxg*nyg*nzg)+'d'),\
            ('qrad_1ray',endian+str(           nxg*nyg*nzg)+'d'),\
            ])
        qqq = np.fromfile(f,dtype=dtype,count=1)
        # MPI position
        ix00 = self.p['xyz'][np0,0]
        jx00 = self.p['xyz'][np0,1]
        kx00 = self.p['xyz'][np0,2]
        
        i0, i1 = nx*ix00, nx*(ix00+1)
        j0, j1 = ny*jx00, ny*(jx00+1)
        k0, k1 = nz*kx00, nz*(kx00+1)
        print(i0,i1,j0,j1,k0,k1)
        
        self.qd['in'][:,:,:,:,i0:i1,j0:j1,k0:k1] = \
            qqq['in'].reshape((mrad,2,2,2,nxg,nyg,nzg),order='F')[:,:,:,:,margin:nxg-margin,margin:nyg-margin,margin:nzg-margin]
        self.qd['qrad'][i0:i1,j0:j1,k0:k1] = \
            qqq['qrad']     .reshape((nxg,nyg,nzg),order='F')[margin:nxg-margin,margin:nyg-margin,margin:nzg-margin]
        self.qd['qrad_1ray'][i0:i1,j0:j1,k0:k1] = \
            qqq['qrad_1ray'].reshape((nxg,nyg,nzg),order='F')[margin:nxg-margin,margin:nyg-margin,margin:nzg-margin]
        f.close()


    # f = open(self.p['datadir']+"qq/qq.dac."+step,'rb')
    # if(order_3D == 1):
    #     self.qc = np.fromfile(f,self.p['endian']+'d',mtype*ixg*jxg*kxg).reshape((mtype,ixg,jxg,kxg),order="F")    
    # if(order_3D == 2):
    #     self.qc = np.fromfile(f,self.p['endian']+'d',mtype*ixg*jxg*kxg).reshape((ixg,mtype,jxg,kxg),order="F")    
    # if(order_3D == 3):
    #     self.qc = np.fromfile(f,self.p['endian']+'d',mtype*ixg*jxg*kxg).reshape((ixg,jxg,mtype,kxg),order="F")    
    # if(order_3D == 4):
    #     self.qc = np.fromfile(f,self.p['endian']+'d',mtype*ixg*jxg*kxg).reshape((ixg,jxg,kxg,mtype),order="F")    

    #f.close()
    
    # if not silent :
    #     print('### variables are stored in self.qc ###')

##############################
def read_qq_slice(self,n_slice,direc,n,silent=False):
    '''
    This method reads 2D data of slice.
    The data is stored in self.ql dictionary

    Parameters:
        n_slice (int): index of slice
        direc (str): slice direction. 'x', 'y', or 'z'
        n (int): a selected time step for data
        silent (bool): True suppresses a message of store
    '''
    import numpy as np

    mtype = self.p['mtype']

    if self.p['geometry'] == 'YinYang':
        files = ['_yin','_yan']
        qls = [self.ql_yin,self.ql_yan]
    else:
        files = ['']
        qls = [self.ql]
    
    for file,ql in zip(files,qls):
        f = open(self.p['datadir']+'slice/qq'+direc+file
                 +'.dac.'+'{0:08d}'.format(n)+'.'
                 +'{0:08d}'.format(n_slice+1),'rb')
        if direc == 'x':
            if self.p['geometry'] == 'YinYang':
                n1, n2 = self.p['jx_yy']+2*self.p['margin'], self.p['kx_yy']+2*self.p['margin']
            else:
                n1, n2 = self.p['jx'], self.p['kx']
        if direc == 'y':
            n1, n2 = self.p['ix'   ], self.p['kx']
        if direc == 'z':
            n1, n2 = self.p['ix'   ], self.p['jx']
        qq_slice = np.fromfile(f,self.p['endian']+'f',(mtype+2)*n1*n2)
    
        ql['ro'] = qq_slice.reshape((n1,n2,mtype+2),order='F')[:,:,0]
        ql['vx'] = qq_slice.reshape((n1,n2,mtype+2),order='F')[:,:,1]
        ql['vy'] = qq_slice.reshape((n1,n2,mtype+2),order='F')[:,:,2]
        ql['vz'] = qq_slice.reshape((n1,n2,mtype+2),order='F')[:,:,3]
        ql['bx'] = qq_slice.reshape((n1,n2,mtype+2),order='F')[:,:,4]
        ql['by'] = qq_slice.reshape((n1,n2,mtype+2),order='F')[:,:,5]
        ql['bz'] = qq_slice.reshape((n1,n2,mtype+2),order='F')[:,:,6]
        ql['se'] = qq_slice.reshape((n1,n2,mtype+2),order='F')[:,:,7]
        ql['ph'] = qq_slice.reshape((n1,n2,mtype+2),order='F')[:,:,8]
        ql['pr'] = qq_slice.reshape((n1,n2,mtype+2),order='F')[:,:,9]
        ql['te'] = qq_slice.reshape((n1,n2,mtype+2),order='F')[:,:,10]
        info = {}
        info['direc'] = direc
        if direc == 'x': slice = self.p['x_slice']
        if direc == 'y': slice = self.p['y_slice']
        if direc == 'z': slice = self.p['z_slice']
        info['slice'] = slice[n_slice]
        info['n_slice'] = n_slice
        ql['info'] = info

        f.close()    
        
    if not silent :
        if self.p['geometry'] == 'YinYang':
            print('### variables are stored in self.ql_yin and self.ql_yan ###')
        else:
            print('### variables are stored in self.ql ###')
    
##############################
def read_qq_2d(self,n,silent=False):
    '''
    This method reads full data of 2D calculation
    The data is stored in self.q2 dictionary

    Parameters:
        n (int): a selected time step for data
        silent (bool): True suppresses a message of store
    '''
    import numpy as np
    
    mtype = self.p["mtype"]
    ix = self.p["ix"]
    jx = self.p["jx"]

    ### Only when memory is not allocated 
    ### and the size of array is different
    ### memory is allocated
    memflag = True
    
    if 'ro' in self.q2:
        memflag = not self.q2['ro'].shape == (ix,jx)
    if 'ro' not in self.q2 or memflag:
        print('memory is newly allocated')
        self.q2["ro"] = np.zeros((ix,jx))
        self.q2["vx"] = np.zeros((ix,jx))
        self.q2["vy"] = np.zeros((ix,jx))
        self.q2["vz"] = np.zeros((ix,jx))
        self.q2["bx"] = np.zeros((ix,jx))
        self.q2["by"] = np.zeros((ix,jx))
        self.q2["bz"] = np.zeros((ix,jx))
        self.q2["se"] = np.zeros((ix,jx))
        self.q2["pr"] = np.zeros((ix,jx))
        self.q2["te"] = np.zeros((ix,jx))
        self.q2["op"] = np.zeros((ix,jx))
        self.q2["tu"] = np.zeros((ix,jx))
        self.q2["fr"] = np.zeros((ix,jx))

    dtyp=np.dtype([ ("qq",self.p["endian"]+str((mtype+5)*ix*jx)+"f") ])
    f = open(self.p['datadir']+"remap/qq/qq.dac."+'{0:08d}'.format(n),'rb')
    qqq = np.fromfile(f,dtype=dtyp,count=1)
    self.q2["ro"] = qqq["qq"].reshape((mtype+5,ix,jx),order="F")[0,:,:]
    self.q2["vx"] = qqq["qq"].reshape((mtype+5,ix,jx),order="F")[1,:,:]
    self.q2["vy"] = qqq["qq"].reshape((mtype+5,ix,jx),order="F")[2,:,:]
    self.q2["vz"] = qqq["qq"].reshape((mtype+5,ix,jx),order="F")[3,:,:]
    self.q2["bx"] = qqq["qq"].reshape((mtype+5,ix,jx),order="F")[4,:,:]
    self.q2["by"] = qqq["qq"].reshape((mtype+5,ix,jx),order="F")[5,:,:]
    self.q2["bz"] = qqq["qq"].reshape((mtype+5,ix,jx),order="F")[6,:,:]
    self.q2["se"] = qqq["qq"].reshape((mtype+5,ix,jx),order="F")[7,:,:]
    self.q2["pr"] = qqq["qq"].reshape((mtype+5,ix,jx),order="F")[8,:,:]
    self.q2["te"] = qqq["qq"].reshape((mtype+5,ix,jx),order="F")[9,:,:]
    self.q2["op"] = qqq["qq"].reshape((mtype+5,ix,jx),order="F")[10,:,:]
    self.q2["tu"] = qqq["qq"].reshape((mtype+5,ix,jx),order="F")[11,:,:]
    self.q2["fr"] = qqq["qq"].reshape((mtype+5,ix,jx),order="F")[12,:,:]
    f.close()

    if not silent :
        print('### variables are stored in self.q2 ###')

def YinYangSet(self):
    '''
    YinYangSet function sets up the YinYang geometry for
    YinYang direct plot. 

    Parameters:
        self (object): R2D2 class instance

    Returns:
        None
    '''
    import numpy as np
    if self.p['geometry'] == 'YinYang':
        if not 'Z_yy' in self.p:
            print('Yes')
            self.p['Z_yy'] , self.p['Y_yy']  = np.meshgrid(self.p['z_yy'],self.p['y_yy'])
            self.p['Zg_yy'], self.p['Yg_yy'] = np.meshgrid(self.p['zg_yy'],self.p['yg_yy'])

            # Geometry in Yang grid        
            self.p['Yo_yy'] = np.arccos(np.sin(self.p['Y_yy'])*np.sin(self.p['Z_yy']))
            self.p['Zo_yy'] = np.arcsin(np.cos(self.p['Y_yy'])/np.sin(self.p['Yo_yy']))
            
            self.p['Yog_yy'] = np.arccos(np.sin(self.p['Yg_yy'])*np.sin(self.p['Zg_yy']))
            self.p['Zog_yy'] = np.arcsin(np.cos(self.p['Yg_yy'])/np.sin(self.p['Yog_yy']))
            
            sct =  np.sin(self.p['Y_yy'])*np.cos(self.p['Z_yy'])
            sco = -np.sin(self.p['Yo_yy'])*np.cos(self.p['Zo_yy'])

            sctg =  np.sin(self.p['Yg_yy'])*np.cos(self.p['Zg_yy'])
            scog = -np.sin(self.p['Yog_yy'])*np.cos(self.p['Zog_yy'])
            
            tmp = sct*sco < 0
            self.p['Zo_yy'][tmp] = np.sign(self.p['Zo_yy'][tmp])*np.pi - self.p['Zo_yy'][tmp]

            tmp = sctg*scog < 0
            self.p['Zog_yy'][tmp] = np.sign(self.p['Zog_yy'][tmp])*np.pi - self.p['Zog_yy'][tmp]

            # self.p['coss_yy'] = -np.sin(self.p['zo_yy'])*np.sin(self.p['zz_yy'])
            # self.p['sins_yy'] =  np.cos(self.p['zo_yy'])/np.sin(self.p['yy_yy'])
            # self.p['cossg_yy'] = -np.sin(self.p['zog_yy'])*np.sin(self.p['zzg_yy'])
            # self.p['sinsg_yy'] =  np.cos(self.p['zog_yy'])/np.sin(self.p['yyg_yy'])
