import numpy as num
import DET_STALTA
import LOC_STALTA

class Stacktraces:

    def __init__(self, tobj, wobj, **inputs):
        # check input objects tobj=traveltime_object wobj=Raw_waveform_object
        self.check_sampling_rate(wobj)
        self.check_starting_time(wobj)

        if ('vfunc' in inputs) or ('hfunc' in inputs):
            vfunc = inputs['vfunc']
            hfunc = inputs['hfunc']
            epsilon = inputs['epsilon']
            derivative = inputs['derivative']
            self.loki_input(wobj, tobj, derivative)
            self.characteristic_function(vfunc, hfunc, epsilon)
            self.enforce_hybrid_ps()
            #self.adaptive_hz_balance()
        else:
            normalize = inputs.get('normthrd', False)
            self.loki_input(wobj, tobj, derivative=False, direct_input=True, normalize=normalize)

            if 'ppower' in inputs:
                self.obs_dataV = self.obs_dataV ** inputs['ppower']
                self.obs_dataH = self.obs_dataH ** inputs['ppower']

    def adaptive_hz_balance(self, Rmax=1.25, eps=1e-12):
        """
        Dynamically rebalance amplitudes if horizontal or vertical energy dominates.
        Acts on characteristic functions (obs_dataH / obs_dataV), using per-sensor RMS.
        """

        if not hasattr(self, "obs_dataH") or not hasattr(self, "obs_dataV"):
            return

        # --- Compute per-station RMS ---
        rmsH = num.sqrt(num.mean(self.obs_dataH**2, axis=1))
        rmsV = num.sqrt(num.mean(self.obs_dataV**2, axis=1))

        # Only include non-zero RMS stations
        validH = rmsH > 0
        validV = rmsV > 0

        if not num.any(validH) or not num.any(validV):
            print("[INFO] Adaptive H/Z balance skipped, only one type of component present")
            return

        mean_rmsH = num.mean(rmsH[validH])
        mean_rmsV = num.mean(rmsV[validV])

        R = mean_rmsH / (mean_rmsV + eps)

        if R > Rmax:
            scaleH = mean_rmsV / (mean_rmsH + eps)
            self.obs_dataH *= scaleH
            print(f"[INFO] Adaptive H/Z balance applied (H dominant, per-station RMS): R={R:.2f}, scaleH={scaleH:.3e}")
        elif 1/R > Rmax:
            scaleV = mean_rmsH / (mean_rmsV + eps)
            self.obs_dataV *= scaleV
            print(f"[INFO] Adaptive H/Z balance applied (V dominant, per-station RMS): 1/R={1/R:.2f}, scaleV={scaleV:.3e}")
        else:
            print(f"[INFO] Adaptive H/Z balance not needed: R={R:.2f}")



    def check_sampling_rate(self, wobj):
        deltas = [wobj.stream[comp][sta][1] for comp in wobj.stream for sta in wobj.stream[comp]]
        deltas = num.array(deltas)
        ideltas = num.unique((deltas * 1E6).astype(int))
        if num.size(ideltas) == 1:
            self.deltat = deltas[0]
        else:
            raise ValueError("Error!! All traces must have the same sampling rate")

    def check_starting_time(self, wobj):
        dtimes = []
        self.ns = 0
        for comp in wobj.stream:
            for sta in wobj.stream[comp]:
                dtimes.append(wobj.stream[comp][sta][0])
                if self.ns < num.size(wobj.stream[comp][sta][2]):
                    self.ns = num.size(wobj.stream[comp][sta][2])
        self.dtime_max = max(dtimes)
        self.evid = self.dtime_max.isoformat()

    def loki_input(self, wobj, tobj, derivative, direct_input=False, normalize=True):
        """
        Load data from waveform object, classify stations, and apply H/Z scaling once.
        """
        if direct_input:
            # Direct assignment for single-component mode
            self.obs_dataV = self.select_data('P', wobj, tobj.db_stations, derivative, normalize)
            self.obs_dataH = self.select_data('S', wobj, tobj.db_stations, derivative, normalize)
            return

        available_comp = list(wobj.stream.keys())
        print(f"[INFO] Available components: {available_comp}")
        self.comp = tuple(available_comp)

        # --- Load raw traces ---
        self.xtr = self.select_data('E', wobj, tobj.db_stations, derivative, normalize)
        self.ytr = self.select_data('N', wobj, tobj.db_stations, derivative, normalize)
        self.ztr = self.select_data('Z', wobj, tobj.db_stations, derivative, normalize)

        # --- Classify stations ---
        self.single_comp_idx = []
        self.three_comp_idx = []
        for i in range(self.nstation):
            comps_present = []
            if num.any(self.xtr[i, :] != 0): comps_present.append('E')
            if num.any(self.ytr[i, :] != 0): comps_present.append('N')
            if num.any(self.ztr[i, :] != 0): comps_present.append('Z')

            if len(comps_present) == 1:
                self.single_comp_idx.append(i)
            elif len(comps_present) == 3:
                self.three_comp_idx.append(i)

        print(f"[INFO] Single-component stations: {self.single_comp_idx}")
        print(f"[INFO] Three-component stations: {self.three_comp_idx}")

        # --- H/Z scaling (per-sensor RMS) ---
        rmsH = num.sqrt(num.mean(self.xtr**2 + self.ytr**2, axis=1))
        rmsZ = num.sqrt(num.mean(self.ztr**2, axis=1))

        validH = rmsH > 0
        validZ = rmsZ > 0

        if num.any(validH) and num.any(validZ):
            mean_rmsH = num.mean(rmsH[validH])
            mean_rmsZ = num.mean(rmsZ[validZ])

            scaleH = mean_rmsZ / (mean_rmsH + 1e-12)
            scaleZ = 1.0

            #self.xtr *= scaleH
            #self.ytr *= scaleH
            #self.ztr *= scaleZ

            self.scaleH = scaleH
            self.scaleZ = scaleZ

            print(f"[INFO] H/Z scaling applied (per-sensor RMS): scaleH={scaleH:.3e}, scaleZ={scaleZ:.3e}")
        else:
            print("[INFO] Skipping H/Z scaling, only one type of component present")
            self.scaleH = 1.0
            self.scaleZ = 1.0

                

    def select_data(self, comp, wobj, db_stations, derivative, normalize):
        self.stations = tuple(wobj.data_stations & set(db_stations))
        self.nstation = num.size(self.stations)
        tr = num.zeros([self.nstation, self.ns])
        stream_comp = wobj.stream.get(comp, {})
        stream_comp = {str(k).strip(): v for k, v in stream_comp.items()}

        for i, sta in enumerate(self.stations):
            if sta not in stream_comp:
                continue
            entry = stream_comp[sta]
            if not (isinstance(entry, (list, tuple)) and len(entry) >= 3):
                continue
            start_dt, delta, data = entry[0], entry[1], entry[2]
            if data is None or len(data) == 0:
                continue
            nstr = num.size(data)
            try:
                idt = int((self.dtime_max - start_dt).total_seconds() / self.deltat)
                if idt < 0: idt = 0
            except Exception:
                idt = 0
            if idt < nstr:
                tr[i, 0:nstr - idt] = data[idt:]
            if derivative:
                tr[i, 1:self.ns] = (tr[i, 1:] - tr[i, :-1]) / self.deltat
                tr[i, 0] = 0.0
            if isinstance(normalize, float):
                trmax = num.max(num.abs(tr[i, :]))
                if trmax >= normalize and trmax > 0:
                    tr[i, :] = tr[i, :] / trmax
            elif normalize:
                trmax = num.max(num.abs(tr[i, :]))
                if trmax > 0:
                    tr[i, :] = tr[i, :] / trmax
        return tr

    # ---------- Analytic signal ----------
    def analytic_signal(self, trace):
        trace = num.atleast_2d(trace)
        nsta, nsamp = trace.shape
        tracef = num.fft.fft(trace, axis=1)
        freqs = num.fft.fftfreq(nsamp, self.deltat)
        traceh = -1j * num.sign(freqs) * tracef
        trace_as = trace + 1j * num.fft.ifft(traceh, axis=1).real
        if trace_as.shape[0] == 1:
            return trace_as[0]
        return trace_as

    # ---------- Characteristic function ----------
    def characteristic_function(self, vfunc='erg', hfunc='pca', epsilon=0.001):
        if len(self.comp) == 1:
            self.cfunc_pca(epsilon)
        else:
            if vfunc == 'erg' and hfunc == 'pca':
                self.cfunc_erg(True)
                self.cfunc_pca(epsilon)
            elif vfunc == 'pca' and hfunc == 'pca':
                self.cfunc_pcafull(epsilon)
            elif vfunc == 'erg' and hfunc == 'erg':
                self.cfunc_erg(False)
            elif vfunc == 'erg' and hfunc == 'null':
                self.cfunc_erg(True)
            elif vfunc == 'cosh' and hfunc == 'cosh':
                self.cfunc_cosh(False)
            elif vfunc == 'tkeo' and hfunc == 'tkeo':
                self.cfunc_tkeo(True)
            else:
                print('wrong characteristic functions, energy used as default')
                self.cfunc_erg(False)

    # ---------- cfunc implementations ----------
    def cfunc_erg(self, ergz):
        if ergz:
            obs_dataV = self.ztr ** 2
            for i in range(self.nstation):
                m = num.max(obs_dataV[i, :])
                if m > 0:
                    obs_dataV[i, :] /= m
            self.obs_dataV = obs_dataV * self.scaleZ
        else:
            obs_dataV = self.ztr ** 2
            obs_dataH = self.xtr ** 2 + self.ytr ** 2
            for i in range(self.nstation):
                mh = num.max(abs(obs_dataH[i, :]))
                mv = num.max(abs(obs_dataV[i, :]))
                if mh > 0:
                    obs_dataH[i, :] /= mh
                if mv > 0:
                    obs_dataV[i, :] /= mv
            self.obs_dataH = obs_dataH * self.scaleH
            self.obs_dataV = obs_dataV * self.scaleZ

    def cfunc_tkeo(self):

        #obs_dataV = self.ztr[1:-1]**2 - self.ztr[:-2]*self.ztr[2:] 
        obs_dataV = self.ztr**2 
        obs_dataH = self.xtr + self.ytr 
        obs_dataH = obs_dataH[1:-1]**2 - obs_dataH[:-2]*obs_dataH[2:]
        for i in range(self.nstation):
            mh = num.max(abs(obs_dataH[i, :]))
            mv = num.max(abs(obs_dataV[i, :]))
            if mh > 0:
                obs_dataH[i, :] /= mh
            if mv > 0:
                obs_dataV[i, :] /= mv

        self.obs_dataH = obs_dataH * self.scaleH
        self.obs_dataV = obs_dataV * self.scaleZ


    def enforce_hybrid_ps(self):
        """
        Enforce P/S usage for hybrid and 1C networks:
        - 3C stations keep true P (V) and S (H)
        - 1C stations use the SAME CF for both P and S
        """

        # If characteristic functions are missing, do nothing
        if not hasattr(self, "obs_dataV") and not hasattr(self, "obs_dataH"):
            return

        # Ensure both arrays exist
        if not hasattr(self, "obs_dataV"):
            self.obs_dataV = num.zeros_like(self.obs_dataH)
        if not hasattr(self, "obs_dataH"):
            self.obs_dataH = num.zeros_like(self.obs_dataV)

        # Loop over single-component stations
        for i in self.single_comp_idx:
            hasV = num.any(self.obs_dataV[i, :] != 0)
            hasH = num.any(self.obs_dataH[i, :] != 0)

            # If only one exists, copy it to the other
            if hasV and not hasH:
                self.obs_dataH[i, :] = self.obs_dataV[i, :].copy()
            elif hasH and not hasV:
                self.obs_dataV[i, :] = self.obs_dataH[i, :].copy()
            # ---- station-count normalization ----
        nz = 0
        nh = 0

        for i in range(self.nstation):
            if num.any(self.obs_dataV[i, :] != 0):
                nz += 1
            if num.any(self.obs_dataH[i, :] != 0):
                nh += 1

        if nz > 0 and nh > 0:
            self.obs_dataV *= (1.0 / nz)
            self.obs_dataH *= (1.0 / nh)


    def cfunc_pca(self, epsilon=0.001):
        obs_dataH = num.zeros([self.nstation, self.ns])

        # Analytic signals
        x_as = self.analytic_signal(self.xtr)
        y_as = self.analytic_signal(self.ytr)
        x_asC = num.conjugate(x_as)
        y_asC = num.conjugate(y_as)

        for i in range(self.nstation):
            for j in range(self.ns):
                cov = num.array([[x_as[i,j]*x_asC[i,j], x_as[i,j]*y_asC[i,j]],
                                [y_as[i,j]*x_asC[i,j], y_as[i,j]*y_asC[i,j]]])
                _, s, _ = num.linalg.svd(cov)
                obs_dataH[i,j] = s[0]**2

            # Per-station normalization for numerical stability
            max_val = num.max(obs_dataH[i, :])
            if max_val > 0:
                obs_dataH[i, :] = (obs_dataH[i, :] / max_val) + epsilon

        # Apply H/Z scaling once at the end
        self.obs_dataH = obs_dataH * self.scaleH


    def cfunc_cosh(self, coshz):
        if coshz:
            obs_dataV = num.cosh(self.ztr)
            for i in range(self.nstation):
                if num.max(obs_dataV[i, :]) > 0:
                    obs_dataV[i, :] /= num.max(obs_dataV[i, :])
            self.obs_dataV = obs_dataV
        else:
            obs_dataV = num.cosh(self.ztr)
            obs_dataH = num.cosh(self.xtr) * num.cosh(self.ytr)
            for i in range(self.nstation):
                if abs(num.max(obs_dataH[i, :])) > 0:
                    obs_dataH[i, :] /= num.max(obs_dataH[i, :])
                if abs(num.max(obs_dataV[i, :])) > 0:
                    obs_dataV[i, :] /= num.max(obs_dataV[i, :])
            self.obs_dataH = obs_dataH * self.scaleH
            self.obs_dataV = obs_dataV * self.scaleZ

    def cfunc_pcafull(self, epsilon=0.001):
        obs_dataH = num.zeros([self.nstation, self.ns])
        obs_dataV = num.zeros([self.nstation, self.ns])

        # Analytic signals
        x_as = self.analytic_signal(self.xtr)
        y_as = self.analytic_signal(self.ytr)
        z_as = self.analytic_signal(self.ztr)
        x_asC = num.conjugate(x_as)
        y_asC = num.conjugate(y_as)
        z_asC = num.conjugate(z_as)

        for i in range(self.nstation):
            for j in range(self.ns):
                # 3D covariance
                cov3d = num.array([
                    [x_as[i,j]*x_asC[i,j], x_as[i,j]*y_asC[i,j], x_as[i,j]*z_asC[i,j]],
                    [y_as[i,j]*x_asC[i,j], y_as[i,j]*y_asC[i,j], y_as[i,j]*z_asC[i,j]],
                    [z_as[i,j]*x_asC[i,j], z_as[i,j]*y_asC[i,j], z_as[i,j]*z_asC[i,j]],
                ])
                # 2D horizontal covariance
                cov2d = cov3d[:2, :2]

                # SVD
                _, s2d, _ = num.linalg.svd(cov2d)
                _, s3d, V3d = num.linalg.svd(cov3d)

                # Horizontal / Vertical decomposition
                obs_dataV[i,j] = (s3d[0] ** 2) * num.abs(V3d[0,2])
                obs_dataH[i,j] = (s2d[0] ** 2) * (1 - num.abs(V3d[0,2]))

            # Per-station normalization for stability
            maxH = num.max(obs_dataH[i, :])
            maxV = num.max(obs_dataV[i, :])
            if maxH > 0:
                obs_dataH[i, :] = (obs_dataH[i, :] / maxH) + epsilon
            if maxV > 0:
                obs_dataV[i, :] /= maxV

        # Apply H/Z scaling once at the end
        self.obs_dataH = obs_dataH * self.scaleH
        self.obs_dataV = obs_dataV * self.scaleZ


    # ---------- loc_stalta ----------
    def loc_stalta(self, nshort_p, nshort_s, slrat, norm=1):
        tshort_p = nshort_p * self.deltat
        tshort_s = nshort_s * self.deltat
        tlong_p = tshort_p * slrat
        tlong_s = tshort_s * slrat
        ks_p = self.deltat / tshort_p
        kl_p = self.deltat / tlong_p
        ks_s = self.deltat / tshort_s
        kl_s = self.deltat / tlong_s
        if len(self.comp) == 1:
            obs_dataP = LOC_STALTA.recursive_stalta(tshort_p, tlong_p, self.deltat, self.obs_dataH, kl_p, ks_p, norm)
            obs_dataS = LOC_STALTA.recursive_stalta(tshort_s, tlong_s, self.deltat, self.obs_dataH, kl_s, ks_s, norm)
        else:
            obs_dataP = LOC_STALTA.recursive_stalta(tshort_p, tlong_p, self.deltat, self.obs_dataV, kl_p, ks_p, norm)
            obs_dataS = LOC_STALTA.recursive_stalta(tshort_s, tlong_s, self.deltat, self.obs_dataH, kl_s, ks_s, norm)
        return obs_dataP * self.scaleZ, obs_dataS * self.scaleH

    # ---------- det_stalta ----------
    def det_stalta(self, nshort_p, nshort_s, slrat, staltap0, staltas0, thres=0.):
        tshort_p = nshort_p * self.deltat
        tshort_s = nshort_s * self.deltat
        tlong_p = tshort_p * slrat
        tlong_s = tshort_s * slrat
        obs_dataP = num.zeros([self.nstation, self.ns])
        obs_dataS = num.zeros([self.nstation, self.ns])
        for i in range(self.nstation):
            obs_dataP[i, :], stltp0 = DET_STALTA.recursive_stalta(staltap0[i][0], staltap0[i][1], tshort_p, tlong_p, self.deltat, self.obs_dataV[i, :], thres)
            obs_dataS[i, :], stlts0 = DET_STALTA.recursive_stalta(staltas0[i][0], staltas0[i][1], tshort_s, tlong_s, self.deltat, self.obs_dataH[i, :], thres)
            staltap0[i][0] = stltp0[0]
            staltap0[i][1] = stltp0[1]
            staltas0[i][0] = stlts0[0]
            staltas0[i][1] = stlts0[1]
        return obs_dataP, obs_dataS, staltap0, staltas0
