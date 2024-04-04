class pj_mapper_v2():
    #-------------------------------------------------------------------#
    def __init__(self,River, parm):
        import matplotlib.cm as cm
        import pandas as pd
        self.river = River
        self.river.mapper = self
        self.parm = parm
        self.colormap = cm.gist_rainbow_r
        self.date1= None; self.date2 = None
        self.ContourLevels = [0,10,10]
        self.vmin = 0; self.vmax = 10
        self.round = 3
        self.min_sample =30
        self.ContourColor = 'grey'
        self.acThreshold= 15 #km2, for network display only
        self.percent = 50; self.stat_func= np.median
        self.map_station_dots = True
        
        self.slideshow_start_text = '1970-1-1'
        self.slideshow_end_text = '2009-1-1'
        self.slideshow_window_text = '1825D'
        self.slideshow_freq_text = '365D'
        self.dpi=72
    #-------------------------------------------------------------------#
        

    def interpolate(self):
        
        interp = pj.interpolator(self.river)
        self.interp = interp
        interp.parm = self.parm
        interp.date1 = self.date1
        interp.date2 = self.date2
        interp.min_sample = self.min_sample
        interp.stat_func = self.stat_func
        interp.kriging()

    #-------------------------------------------------------------------#
        
    def plot_network(self):
        
        River = self.river
        threshold = self.acThreshold
        rt = River.nodes[0]
        
        # plot matrix #
        plt.matshow(self.interp.z_matrix,cmap=self.colormap, origin='lower', 
                    aspect='auto',vmin=self.vmin, vmax=self.vmax)
        plt.gca().xaxis.tick_bottom()
        plt.colorbar()
        
        # plot sampled stations #
        if self.map_station_dots:
            for st in River.stations:
                if not(st.sampled):
                    continue
                if st.source == 'PWQMN':
                    plt.scatter(np.sqrt(st.area()),st.stream.pjNum2(),
                                marker="o",edgecolors='black',
                                color="orange",s=50, zorder=4, alpha = 1, 
                                edgecolor = 'black' )
                    
                elif st.source == "LONDON":
                    plt.scatter(np.sqrt(st.area()),st.stream.pjNum2(),
                                marker="o",edgecolors='black',color="red",
                                s=50, zorder=4, alpha = 1, edgecolor = 'black' )
        
        # Plot stream lines #

        for i in rt.iter_pj():
            if i.ac < threshold :
                continue

            if i.rc:
                if i.rc.ac> threshold :
                    lineR = line([np.sqrt(i.ac),i.pjNum2()],[np.sqrt(i.rc.ac),i.rc.pjNum2()])

                    plt.plot([lineR.x1,lineR.x2],[lineR.y1,lineR.y2],color='w', zorder=1)

            if i.lc:
                if i.lc.ac> threshold :
                    tailAc =River.flow_acc[i.cells[-1][0],i.cells[-1][1]]
                    acRatio =1- River.flow_acc[i.cells[-1][0],i.cells[-1][1]]/float(i.ac)
                    fork = lineR.position(acRatio)
                    lineL = line(fork,[np.sqrt(i.lc.ac),i.lc.pjNum2()])
                    plt.plot([lineL.x1,lineL.x2],[lineL.y1,lineL.y2],color = 'w', zorder=1)

        # final adjustments #
        plt.gca().tick_params(axis = 'both', which = 'major', labelsize = 20)
        plt.ylim(ymin=0)
        plt.gcf().set_size_inches(10,12)
        
        
        self.fig = plt.gcf()
        
    #-------------------------------------------------------------------#
        
    def trim(self):
        from scipy import interpolate
        
        cols = self.interp.z_matrix.shape[1]
        rows = self.interp.z_matrix.shape[0]
        
        main_channel = self.river.nodes[0].iter_pj()[-1].downstreams()
        xs = np.array([nd.ac**0.5 for nd in main_channel])
        ys = np.array([nd.pjNum2() for nd in main_channel])
        f = interpolate.interp1d(xs, ys+5, fill_value='extrapolate',kind='cubic')

        xs2 = np.arange(0,max(self.interp.xs)+20,0.1)
        ys2 = f(xs2)
        
        idx = len(ys2[ys2>=0])
        ys2 = ys2[:idx]
        xs2 = xs2[:idx]

        xs2 = np.append(xs2,[xs2[-1], cols, cols, -0.5, -0.5,   xs2[0]])
        ys2 = np.append(ys2,[-0.5,    -0.5, rows, rows, ys2[0], ys2[0]])
        
        plt.fill(xs2,ys2,color='#FBFBFB',zorder=2)
        plt.gca().margins(x=0)
        plt.gca().margins(y=0)

        plt.tight_layout()
        self.fig = plt.gcf()
        self.ax = plt.gca()
    #-----------------------------------------------------------------------#
    def slideshow_estimate_max_min(self):
        min_sample = self.min_sample
        self.slideshow_window = pd.to_timedelta(self.slideshow_window_text)
        stat_func= self.stat_func
        df = self.river[self.parm]
        
        mean_list = []
        for dt in pd.date_range(start=self.slideshow_start_text, 
                                end= self.slideshow_end_text, 
                                freq=self.slideshow_freq_text):
            date1 = dt
            date2 = date1 + self.slideshow_window
            mask = (df['DATE'] > date1) & (df['DATE'] < date2)

            df2 = df[mask]
            stns = df2['STATION'].value_counts()[df2['STATION'].value_counts()>min_sample].index
            df3 = df2[df2['STATION'].isin(stns)]

            mean_value = df3.groupby('STATION')['RESULT'].apply(stat_func).tolist()
            mean_list += mean_value

        self.vmin = pd.Series(mean_list).quantile(0.05)
        self.vmax =pd.Series(mean_list).quantile(0.85)
        print self.vmin, self.vmax
    #------------------------------------------------------------------------##
    def locate_row_col(self, x, y, fig, ax):
        
        x, y = ax.transData.transform([[x],[y]]).flatten()
        width, height = fig.canvas.get_width_height()
        row = height-y
        col = x
        return  600.0*row/height, 500.0*col/width
     #------------------------------------------------------------------------##
    
    def plot_slideshow(self):
        import random
        import string
        import os   
        # setup workspace
        uniname = self.parm + '-' + ''.join(random.choice('ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789')
                                            for x in range(6))
        os.makedirs("exports/fig/slideshow/{}".format(uniname))
        os.makedirs("exports/fig/slideshow/{}/slides".format(uniname))
        self.slideshow_window = pd.to_timedelta(self.slideshow_window_text)
        
        
        nIndex = 0
        list_date = []
        list_dqi = []
        frames_obj = {}
        print "generate map series...."
        for dt in pd.date_range(start=self.slideshow_start_text, 
                                end= self.slideshow_end_text, 
                                freq=self.slideshow_freq_text):
            
            try:
                self.date1 = dt
                self.date2 = dt + self.slideshow_window
                self.interpolate()
                self.plot_network()
                self.trim()
                self.fig = plt.gcf()
                self.ax = plt.gca()
                self.fig.dpi = self.dpi
                
                plt.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.1)
                
                self.fig.savefig("exports/fig/slideshow/{}/slides/{}.png".format(
                    uniname, nIndex ), dpi=self.dpi)
                print "w-h: " , self.fig.canvas.get_width_height()
                
                plt.close()
                
                ## log data for json file
                list_date.append((dt+self.slideshow_window/2.0).strftime('%Y-%m-%d'))
                list_dqi.append(np.median(self.interp.cts))
                
                fr_obj = {}
                fr_obj['index'] = nIndex
                fr_obj['active_station_IDs'] =[stn.ID for stn in  self.river.stations if stn.sampled]
                fr_obj['DQI'] = np.median(self.interp.cts)
                fr_obj['date'] = (dt+self.slideshow_window/2.0).strftime('%Y-%m-%d')
                frames_obj[nIndex] = fr_obj  
                nIndex  += 1
                
            except:
                print 'error',dt
                pass

        ## Write index.json file
        print "wirting index.json file"
        import json
        sample = {'id':uniname , 'river': self.river.name, 
                  'parm':self.parm, 'dates':list_date,
                  "numFrame": nIndex, 'DQI':list_dqi
                 }
        
        stations_obj = {}
        for stn in self.river.stations:
            stn_obj = {}
            stn_obj['lat'] = stn.lat_lon[0]
            stn_obj['lon'] = stn.lat_lon[1]
            stn_obj['stream_name'] = stn.name[0]
            stn_obj['location_name'] = stn.name[1]
            stn_obj['ID'] = stn.ID
            stn_obj['area'] = stn.area()
            stn_obj['x'] = np.sqrt(stn.area())
            stn_obj['y'] = stn.stream.pjNum2()
            stations_obj[stn.ID] = stn_obj
            row, col = self.locate_row_col(np.sqrt(stn.area()),
                                      stn.stream.pjNum2(),
                                      self.fig, self.ax)
            #print "writing w-h: " , self.fig.canvas.get_width_height()
            
            stn_obj['row'] = row
            stn_obj['col'] = col
                
            # ------- #
            
        sample['stations']  = stations_obj
        sample['frames']  = frames_obj

        with open("exports/fig/slideshow/{}/info.json".format(uniname), 'w') as fp:
            text = 'data = ' + json.dumps(obj=sample,indent=4)
            fp.write(text)
            