import numpy as num
import DET_STALTA     # C
import LOC_STALTA     # C


class Stacktraces:

    def __init__(self, tobj, wobj, **inputs):
        #check input objects tobj=traveltime_object wobj=Raw_waveform_object
        self.check_sampling_rate(wobj)
        self.check_starting_time(wobj)
        
        if ('vfunc' in inputs) or ('hfunc' in inputs):
            # need to calculate characteristic functions from input data
            vfunc=inputs['vfunc']
            hfunc=inputs['hfunc']
            epsilon=inputs['epsilon']
            derivative=inputs['derivative']
            self.loki_input(wobj, tobj, derivative)
            self.characteristic_function(vfunc, hfunc, epsilon)
        else:
            # no need to calculate characteristic function from input data
            # directly input data of P and S component
            if 'normthrd' in inputs:
                normalize = inputs['normthrd']
            else:
                normalize = False
            self.loki_input(wobj, tobj, derivative=False, direct_input=True, normalize=normalize)
            
            # compute array element wise power over the input probabilities if needed
            if 'ppower' in inputs:
                self.obs_dataV = self.obs_dataV**inputs['ppower']
                self.obs_dataH = self.obs_dataH**inputs['ppower']


    def check_sampling_rate(self,wobj):
        intsamp=1E6
        deltas=[]
        #print(wobj.stream)
        for comp in (wobj.stream).keys():
            #print('comp',comp)
            for sta in (wobj.stream[comp]).keys():
                deltas.append(wobj.stream[comp][sta][1])
        #print('deltas', deltas)
        deltas=num.array(deltas)
        ideltas=num.unique((deltas*intsamp).astype(int))
        if num.size(ideltas)==1:
            self.deltat=deltas[0]
        else:
            raise ValueError('Error!! All trace must have the same sampling rate')


    def check_starting_time(self,wobj):
        dtimes=[]
        self.ns=0
        for comp in (wobj.stream).keys():
            for sta in (wobj.stream[comp]).keys():
                dtimes.append(wobj.stream[comp][sta][0])
                if self.ns<num.size(wobj.stream[comp][sta][2]):
                    self.ns=num.size(wobj.stream[comp][sta][2])
        self.dtime_max=max(dtimes)
        self.evid=(self.dtime_max).isoformat()


    def loki_input(self, wobj, tobj, derivative, direct_input=False, normalize=True):
            if direct_input:
                # directly use input data as characteristic function
                self.obs_dataV = self.select_data('P', wobj, tobj.db_stations, derivative, normalize)
                self.obs_dataH = self.select_data('S', wobj, tobj.db_stations, derivative, normalize)
                return

            # Build xtr, ytr, ztr arrays from whichever components are available.
            # wobj.stream has keys like 'E','N','Z' (each mapping station -> [dtime, delta, data])
            available_comp = list(wobj.stream.keys())
            self.comp = tuple(available_comp)

            # fetch arrays for each component if present, else zeros
            # select_data returns array shape (nstation, ns)
            # We'll use names xtr -> E, ytr -> N, ztr -> Z
            # If a component is missing, select_data will return zeros.
            self.xtr = self.select_data('E', wobj, tobj.db_stations, derivative, normalize)
            self.ytr = self.select_data('N', wobj, tobj.db_stations, derivative, normalize)
            self.ztr = self.select_data('Z', wobj, tobj.db_stations, derivative, normalize)

            # If only one component exists overall, follow original behavior: make ytr same as xtr
            # (this retains older code semantics where single-component used as both x and y for horizontal)
            if len(available_comp) == 1:
                # find which comp is present and copy to both xtr and ytr if needed
                only = available_comp[0]
                if only == 'E':
                    self.xtr = self.select_data('E', wobj, tobj.db_stations, derivative, normalize)
                    self.ytr = self.xtr.copy()
                    self.ztr = num.zeros_like(self.xtr)
                elif only == 'N':
                    self.ytr = self.select_data('N', wobj, tobj.db_stations, derivative, normalize)
                    self.xtr = self.ytr.copy()
                    self.ztr = num.zeros_like(self.ytr)
                elif only == 'Z':
                    self.ztr = self.select_data('Z', wobj, tobj.db_stations, derivative, normalize)
                    # single vertical: set horizontals to zeros
                    self.xtr = num.zeros_like(self.ztr)
                    self.ytr = num.zeros_like(self.ztr)

    def select_data(self, comp, wobj, db_stations, derivative, normalize):
        """
        More robust selection:
        - comp: component letter e.g. 'E','N','Z'
        - if that comp is missing for a station, its row stays zeros
        - aligns traces in time using dtime_max and deltat
        """
        self.stations = tuple(wobj.data_stations & set(db_stations))
        self.nstation = num.size(self.stations)
        tr = num.zeros([self.nstation, self.ns])
        # get dict for this component, or empty dict if comp not present
        stream_comp = wobj.stream.get(comp, {})

        # normalize keys
        stream_comp = {str(k).strip(): v for k, v in stream_comp.items()}

        for i, sta in enumerate(self.stations):
            if sta not in stream_comp:
                # leave zeros (component missing for this station)
                # print only if you want verbose info:
                # print(f"[INFO] Station {sta} missing component {comp}, filling zeros.")
                continue

            entry = stream_comp[sta]
            if not (isinstance(entry, (list, tuple)) and len(entry) >= 3):
                print(f"[WARN] Unexpected entry format for station {sta} comp {comp}: {entry}")
                continue

            start_dt, delta, data = entry[0], entry[1], entry[2]
            if data is None or len(data) == 0:
                print(f"[WARN] Station {sta} component {comp} has empty data, skipping.")
                continue

            nstr = num.size(data)
            # time alignment: number of samples to skip on left
            try:
                idt = int((self.dtime_max - start_dt).total_seconds() / self.deltat)
                if idt < 0:
                    # if this trace starts after reference, clip negative index
                    idt = 0
            except Exception as e:
                print(f"[WARN] Time alignment failed for {sta} comp {comp}: {e}")
                idt = 0

            if idt < nstr:
                tr[i, 0: nstr - idt] = data[idt:]
            # else remain zeros (trace fully before dtime_max)

            if derivative:
                tr[i, 1:self.ns] = (tr[i, 1:] - tr[i, :-1]) / self.deltat
                tr[i, 0] = 0.0

            # normalization options
            if isinstance(normalize, float):
                trmax = num.max(num.abs(tr[i, :]))
                if trmax >= normalize and trmax > 0:
                    tr[i, :] = tr[i, :] / trmax
            elif normalize:
                trmax = num.max(num.abs(tr[i, :]))
                if trmax > 0:
                    tr[i, :] = tr[i, :] / trmax

        return tr

    def analytic_signal(self, trace):
        tracef=num.fft.fft(trace)
        nsta,nfreq=num.shape(tracef)
        freqs=num.fft.fftfreq(nfreq,self.deltat)
        traceh=-1j*num.sign(freqs).T*tracef
        trace=trace+1j*num.fft.ifft(traceh).real
        return trace

    def time_extractor(self, tp, ts):
        nxyz= num.size(tp[self.stations[0]])
        tp_mod=num.zeros([nxyz,self.nstation])
        ts_mod=num.zeros([nxyz,self.nstation])
        for i,sta in enumerate(self.stations):
            tp_mod[:,i]=tp[sta]
            ts_mod[:,i]=ts[sta]
        return (tp_mod, ts_mod)

    def characteristic_function(self, vfunc='erg', hfunc='pca', epsilon=0.001):
        if len(self.comp)==1: 
            #self.cfunc_erg(True)
            self.cfunc_pca(epsilon)

        else:
            if vfunc=='erg' and hfunc=='pca':
                self.cfunc_erg(True)
                self.cfunc_pca(epsilon)
            elif vfunc=='pca' and hfunc=='pca':
                self.cfunc_pcafull(epsilon)
            elif vfunc=='erg' and hfunc=='erg':
                self.cfunc_erg(False)
            elif vfunc=='erg' and hfunc=='null':
                self.cfunc_erg(True)
            elif vfunc=='cosh' and hfunc=='cosh':
                self.cfunc_cosh(False)
            else:
                print('wrong characterstic functions, energy used as default')
                self.cfunc_erg(False)

    def cfunc_erg(self, ergz):
        # ergz True : use only vertical (ztr) as obs_dataV
        # ergz False: compute vertical energy and horizontal energy from available horizontals
        if ergz:
            obs_dataV = (self.ztr ** 2)
            for i in range(self.nstation):
                m = num.max(obs_dataV[i, :])
                if m > 0:
                    obs_dataV[i, :] = obs_dataV[i, :] / m
            self.obs_dataV = obs_dataV
        else:
            # vertical
            obs_dataV = (self.ztr ** 2)
            # horizontal: sum squares of available horizontal components
            # if only E or only N present, the missing one is zeros (per loki_input)
            obs_dataH = (self.xtr ** 2) + (self.ytr ** 2)
            for i in range(self.nstation):
                mh = num.max(abs(obs_dataH[i, :]))
                mv = num.max(abs(obs_dataV[i, :]))
                if mh > 0:
                    obs_dataH[i, :] = obs_dataH[i, :] / mh
                if mv > 0:
                    obs_dataV[i, :] = obs_dataV[i, :] / mv
            self.obs_dataH = obs_dataH
            self.obs_dataV = obs_dataV

    def cfunc_pca(self, epsilon):
        # PCA on horizontal pair (xtr,ytr) — if one horizontal is zeros, PCA will gracefully reduce
        obs_dataH = num.zeros([self.nstation, self.ns])
        obs_dataH1 = self.analytic_signal(self.xtr)
        obs_dataH2 = self.analytic_signal(self.ytr)
        obs_dataH1C = num.conjugate(obs_dataH1)
        obs_dataH2C = num.conjugate(obs_dataH2)
        xx = obs_dataH1 * obs_dataH1C
        xy = obs_dataH1 * obs_dataH2C
        yx = obs_dataH2 * obs_dataH1C
        yy = obs_dataH2 * obs_dataH2C
        for i in range(self.nstation):
            for j in range(self.ns):
                cov = num.array([[xx[i, j], xy[i, j]], [yx[i, j], yy[i, j]]])
                U, s, V = num.linalg.svd(cov, full_matrices=True)
                obs_dataH[i, j] = (s[0] ** 2)
            if abs(num.max(obs_dataH[i, :])) > 0:
                obs_dataH[i, :] = (obs_dataH[i, :] / num.max(obs_dataH[i, :])) + epsilon
        self.obs_dataH = obs_dataH

    def cfunc_cosh(self, coshz):

        if coshz:
            obs_dataV=(num.cosh(self.ztr))
            for i in range(self.nstation):
                if num.max(obs_dataV[i,:]) > 0:
                    obs_dataV[i,:]=(obs_dataV[i,:]/num.max(obs_dataV[i,:]))
            self.obs_dataV=obs_dataV
        else:
            obs_dataV=num.cosh(self.ztr)
            obs_dataH=num.cosh(self.xtr)*num.cosh(self.ytr)
            for i in range(self.nstation):
                if abs(num.max(obs_dataH[i,:])) > 0:
                    obs_dataH[i,:]=(obs_dataH[i,:]/num.max(obs_dataH[i,:]))
                if abs(num.max(obs_dataV[i,:])) > 0:
                    obs_dataV[i,:]=(obs_dataV[i,:]/num.max(obs_dataV[i,:]))
            self.obs_dataH=obs_dataH
            self.obs_dataV=obs_dataV

    def cfunc_pcafull(self, epsilon):
        obs_dataH=num.zeros([self.nstation,self.ns]); obs_dataV=num.zeros([self.nstation,self.ns])
        obs_dataH1=self.analytic_signal(self.xtr); obs_dataH2=self.analytic_signal(self.ytr); obs_dataH3=self.analytic_signal(self.ztr)
        obs_dataH1C=num.conjugate(obs_dataH1); obs_dataH2C=num.conjugate(obs_dataH2); obs_dataH3C=num.conjugate(obs_dataH3)
        xx=obs_dataH1*obs_dataH1C; xy=obs_dataH1*obs_dataH2C; xz=obs_dataH1*obs_dataH3C
        yx=obs_dataH2*obs_dataH1C; yy=obs_dataH2*obs_dataH2C; yz=obs_dataH2*obs_dataH2C
        zx=obs_dataH3*obs_dataH1C; zy=obs_dataH3*obs_dataH2C; zz=obs_dataH3*obs_dataH3C
        for i in range(self.nstation):
            for j in range(self.ns):
                cov3d=num.array([[xx[i,j], xy[i,j], xz[i,j]],[yx[i,j], yy[i,j], yz[i,j]],[zx[i,j],zy[i,j],zz[i,j]]])
                cov2d=num.array([[xx[i,j], xy[i,j]],[yx[i,j], yy[i,j]]])
                U2d, s2d, V2d = num.linalg.svd(cov2d, full_matrices=True)
                U3d, s3d, V3d = num.linalg.svd(cov3d, full_matrices=True)
                obs_dataV[i,j]=(s3d[0]**2)*(num.abs(V3d[0][2]))
                obs_dataH[i,j]=(s2d[0]**2)*(1-num.abs(V3d[0][2]))
                
            if abs(num.max(obs_dataH[i,:])) > 0:
                obs_dataH[i,:]=(obs_dataH[i,:]/num.max(obs_dataH[i,:]))+epsilon
            if abs(num.max(obs_dataV[i,:])) > 0:
                obs_dataV[i,:]=(obs_dataV[i,:]/num.max(obs_dataV[i,:]))
        self.obs_dataH=obs_dataH
        self.obs_dataV=obs_dataV



    # def pca_win(xtr, ytr, ztr, iwin):
    #     for i in range(self.nstation):
    #          for j in range(iwin,self.ns-iwin):
    #              X=self.xtr[i,j-iwin:j+iwin]-num.mean(self.xtr[i,j-iwin:j+iwin])
    #              Y=self.ytr[i,j-iwin:j+iwin]-num.mean(self.ytr[i,j-iwin:j+iwin])
    #              Z=self.ztr[i,j-iwin:j+iwin]-num.mean(self.ztr[i,j-iwin:j+iwin])
    #              cov=num.vstack((X,Y,Z))
    #              C=num.dot(cov,cov.T)
    #              U, s, V = num.linalg.svd(C, full_matrices=True)
    #              obs_dataH[i,j]=1-((s[1]+s[2])/(2*s[0])))
    #              azm.append(num.arctan2(V[0][1],V[0][0])*(180/num.pi))
    #              inc.append(num.arccos(num.abs(V[0][2]))*(180/num.pi))
    #           obs_dataH[i,j]
    #     dol=num.array(dol); azm=num.array(azm); inc=num.array(inc)
    #     azm[azm<0]=360+azm[azm<0]
    # return dol,azm,inc


    def loc_stalta(self, nshort_p, nshort_s, slrat, norm=1):
        tshort_p=nshort_p*self.deltat; tshort_s=nshort_s*self.deltat
        tlong_p=tshort_p*slrat; tlong_s=tshort_s*slrat
        ks_p=self.deltat/tshort_p; kl_p=self.deltat/tlong_p;
        ks_s=self.deltat/tshort_s; kl_s=self.deltat/tlong_s;
        if len(self.comp)==1:
            obs_dataP=LOC_STALTA.recursive_stalta(tshort_p, tlong_p, self.deltat, self.obs_dataH, kl_p, ks_p, norm)
            obs_dataS=LOC_STALTA.recursive_stalta(tshort_s, tlong_s, self.deltat, self.obs_dataH, kl_s, ks_s, norm)
        else: 
            obs_dataP=LOC_STALTA.recursive_stalta(tshort_p, tlong_p, self.deltat, self.obs_dataV, kl_p, ks_p, norm)
            obs_dataS=LOC_STALTA.recursive_stalta(tshort_s, tlong_s, self.deltat, self.obs_dataH, kl_s, ks_s, norm)
        return obs_dataP, obs_dataS


    def det_stalta(self, nshort_p, nshort_s, slrat, staltap0, staltas0, thres=0.):
        tshort_p=nshort_p*self.deltat; tshort_s=nshort_s*self.deltat
        tlong_p=tshort_p*slrat; tlong_s=tshort_s*slrat
        obs_dataP=num.zeros([self.nstation,self.ns])
        obs_dataS=num.zeros([self.nstation,self.ns])
        for i in range(self.nstation):
            obs_dataP[i,:],stltp0=DET_STALTA.recursive_stalta(staltap0[i][0], staltap0[i][1], tshort_p, tlong_p, self.deltat, self.obs_dataV[i,:], thres)
            obs_dataS[i,:],stlts0=DET_STALTA.recursive_stalta(staltas0[i][0], staltas0[i][1], tshort_s, tlong_s, self.deltat, self.obs_dataH[i,:], thres)
            staltap0[i][0]=stltp0[0]; staltap0[i][1]=stltp0[1]
            staltas0[i][0]=stlts0[0]; staltas0[i][1]=stlts0[1]
        return obs_dataP, obs_dataS, staltap0, staltas0