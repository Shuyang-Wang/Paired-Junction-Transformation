from osgeo import gdal
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
from osgeo import osr


## Define global constant

dX = [1, 1, 1, 0, -1, -1, -1, 0]
dY = [-1, 0, 1, 1, 1, 0, -1, -1]
dYX = list(zip(dY,dX))
D8FlowVals = [128, 1, 2, 4, 8, 16, 32, 64]
inFlowVals = [8, 16, 32, 64, 128, 1, 2, 4]
resolution = 30

river_latlon = {
    'Credit':(-79.5833773995,43.5493864802),
    'Thames':(-82.45372694, 42.31865166),
    'Grand':(-79.576616, 42.857583),
    'Saugeen':(-81.356047, 44.502847),
    'Medway':(-81.270671,43.012422),
    'Humber':(-79.475778, 43.634627),
}



def rc_to_latlon(River,r_c_list):
    """
    Convert a list of [row, col] to [lat, lon]
    """
    ds = River.ds
    gt = ds.GetGeoTransform()
    output = []

    old_cs = osr.SpatialReference()
    old_cs.ImportFromWkt(ds.GetProjectionRef())
    wgs84_wkt = """
    GEOGCS["WGS 84",
        DATUM["WGS_1984",
            SPHEROID["WGS 84",6378137,298.257223563,
                AUTHORITY["EPSG","7030"]],
            AUTHORITY["EPSG","6326"]],
        PRIMEM["Greenwich",0,
            AUTHORITY["EPSG","8901"]],
        UNIT["degree",0.01745329251994328,
            AUTHORITY["EPSG","9122"]],
        AUTHORITY["EPSG","4326"]]"""
    new_cs = osr.SpatialReference()
    new_cs.ImportFromWkt(wgs84_wkt)
    transform = osr.CoordinateTransformation(old_cs,new_cs)
    
    for r_c in r_c_list:
        row =River.extent[0] + r_c[0]
        col =River.extent[2] + r_c[1]
        y = gt[3] - row*30
        x = gt[0] + col*30
        LATLON = transform.TransformPoint(x, y)[::-1][1:]
        LATLON = (round(LATLON[0],4), round(LATLON[1],4)) # round to reduce kml size
        output.append(LATLON)
        
    return output
    ## END OF rc_to_latlon


def Flow_Acc(input_pointer):

    '''
    Globel Constant: D8FlowVa
    '
    ls,
    '''

    rows = input_pointer.shape[0]
    cols = input_pointer.shape[1]
    
    def inflow_count(input_pointer):
        outputGrid = np.zeros((rows, cols),dtype=float)
        for i in D8FlowVals:
            shift = dYX[D8FlowVals.index(i)]
            tmpGrid= input_pointer == i
            tmpGrid = tmpGrid.astype(int)
            tmpGrid = np.roll(tmpGrid,shift[0], 0)
            tmpGrid = np.roll(tmpGrid,shift[1], 1)

            if shift[0] == 1:
                tmpGrid[0] = np.linspace(0, 0, len(tmpGrid[0]))
            elif shift[0] == -1:
                tmpGrid[-1] = np.linspace(0, 0, len(tmpGrid[0]))
            if shift[1] == 1:
                tmpGrid[:, 0] = np.linspace(0, 0, tmpGrid.shape[0])
            elif shift[1] == -1:
                tmpGrid[:, -1] = np.linspace(0, 0, tmpGrid.shape[0])
            outputGrid = outputGrid + tmpGrid
        outputGrid[input_pointer<0] =-1
        
        return outputGrid
    
    
    
    flow_accum_raster = np.full((rows, cols), 1, dtype=int)
    tmpGrid = inflow_count(input_pointer)

    for row in range(rows):
        for col in range(cols):
            if tmpGrid[row, col] == 0:
                
                tmpGrid[row, col] = -1
                flag = True
                x = col
                y = row
                while flag:
                    
                    z = flow_accum_raster[y, x]
                    flowDir = input_pointer[y, x]
                    if flowDir > 0:
                        i = D8FlowVals.index(flowDir)
                        x += dX[i]
                        y += dY[i]
                        if y < rows and x < cols:
                            flow_accum_raster[y, x] += z
                            tmpGrid[y, x] -= 1
                            numInFlows = tmpGrid[y, x]
                            if numInFlows == 0:
                                flag = True
                                tmpGrid[y, x] = -1
                            else:
                                flag = False
                        else:
                            flag = False
                    else:
                        flag = False
                        ###
                    ###
                ###
            ###
        ###
    flow_accum_raster = flow_accum_raster * float(resolution)**2/(1000000.0)
    return flow_accum_raster
    ## END OF Flow_Acc()
##
def resi_time(input_pointer,flow_acc):

    '''
    generate average residence time raster
    input: flow direction raster, ArcGIS style
    requirment: self.flow_acc(flow accumulation) has been assigned
    '''

    rows = input_pointer.shape[0]
    cols = input_pointer.shape[1]
    def inflow_count(input_pointer):
        '''
        generate count raster (number of neigbouring infows per cell)
        '''
        outputGrid = np.zeros((rows, cols),dtype=float)
        for i in D8FlowVals:
            shift = dYX[D8FlowVals.index(i)]
            tmpGrid= input_pointer == i
            tmpGrid = tmpGrid.astype(int)
            tmpGrid = np.roll(tmpGrid,shift[0], 0)
            tmpGrid = np.roll(tmpGrid,shift[1], 1)

            if shift[0] == 1:
                tmpGrid[0] = np.linspace(0, 0, len(tmpGrid[0]))
            elif shift[0] == -1:
                tmpGrid[-1] = np.linspace(0, 0, len(tmpGrid[0]))
            if shift[1] == 1:
                tmpGrid[:, 0] = np.linspace(0, 0, tmpGrid.shape[0])
            elif shift[1] == -1:
                tmpGrid[:, -1] = np.linspace(0, 0, tmpGrid.shape[0])

            outputGrid = outputGrid + tmpGrid
        outputGrid[input_pointer<0] =-1
        
        return outputGrid
    
    # initiate, define raster grids
    resi_time_raster = np.full((rows, cols), 1, dtype=float) 
    # residence time, total
    count_raster = np.array(np.round(flow_acc/(resolution**2)*1000000)) 
    # total counts of upstream cells
    tmpGrid = inflow_count(input_pointer) 
    # numuber of inflow for each cell

    for row in range(rows):
        for col in range(cols):

            # if (row % 100 ==0)&(col==1):
            #     print row

            if tmpGrid[row, col] == 0:
                tmpGrid[row, col] = -1
                flag = True
                x = col
                y = row

                while flag:
                    z = resi_time_raster[y, x]
                    flowDir = input_pointer[y, x]
                    
                    if flowDir > 0:
                        i = D8FlowVals.index(flowDir)
                        count = count_raster[y,x]
                        x += dX[i]
                        y += dY[i]

                        if i in [2,8,32,128]:
                            dis = 1.414213562
                        else:
                            dis = 1
                        if y < rows and x < cols:
                            resi_time_raster[y, x] += (z + count*dis)
                            tmpGrid[y, x] -= 1
                            
                            if tmpGrid[y, x] == 0:
                                flag = True
                                tmpGrid[y, x] = -1
                            else:
                                flag = False
                        else:
                            flag = False
                    else:
                        flag = False
                        ###
                    ###
                ###
            ###
        ###
    resi_time_raster = resi_time_raster/count_raster
    return resi_time_raster
    ## END OF resi_time()


def Pop_Acc(input_pointer,pop_density):

    '''
    Globel Constant: D8FlowVa
    '
    ls,
    '''

    rows = input_pointer.shape[0]
    cols = input_pointer.shape[1]
    
    def inflow_count(input_pointer):
        outputGrid = np.zeros((rows, cols),dtype=float)
        for i in D8FlowVals:
            shift = dYX[D8FlowVals.index(i)]
            tmpGrid= input_pointer == i
            tmpGrid = tmpGrid.astype(int)
            tmpGrid = np.roll(tmpGrid,shift[0], 0)
            tmpGrid = np.roll(tmpGrid,shift[1], 1)

            if shift[0] == 1:
                tmpGrid[0] = np.linspace(0, 0, len(tmpGrid[0]))
            elif shift[0] == -1:
                tmpGrid[-1] = np.linspace(0, 0, len(tmpGrid[0]))
            if shift[1] == 1:
                tmpGrid[:, 0] = np.linspace(0, 0, tmpGrid.shape[0])
            elif shift[1] == -1:
                tmpGrid[:, -1] = np.linspace(0, 0, tmpGrid.shape[0])
            outputGrid = outputGrid + tmpGrid
        outputGrid[input_pointer<0] =-1
        
        return outputGrid
    

    pop_accum_raster = pop_density.copy()
    tmpGrid = inflow_count(input_pointer)

    for row in range(rows):
        for col in range(cols):
            if tmpGrid[row, col] == 0:
                
                tmpGrid[row, col] = -1
                flag = True
                x = col
                y = row
                while flag:
                    
                    z = pop_accum_raster[y, x]
                    flowDir = input_pointer[y, x]
                    if flowDir > 0:
                        i = D8FlowVals.index(flowDir)
                        x += dX[i]
                        y += dY[i]
                        if y < rows and x < cols:
                            pop_accum_raster[y, x] += z
                            tmpGrid[y, x] -= 1
                            numInFlows = tmpGrid[y, x]
                            if numInFlows == 0:
                                flag = True
                                tmpGrid[y, x] = -1
                            else:
                                flag = False
                        else:
                            flag = False
                    else:
                        flag = False
                        ###
                    ###
                ###
            ###
        ###
    pop_accum_raster = pop_accum_raster
    return pop_accum_raster
    ## END OF Flow_Acc()
##

def read_raster(fname,dtype):
    
    tmpDS = gdal.Open(fname)
    rows = tmpDS.RasterYSize 
    cols = tmpDS.RasterXSize
    output = tmpDS.ReadAsArray(0, 0, cols, rows).astype(dtype)
    return output
    ## END OF read_raster()
##
def left_top_corner(ds,extent):

    gt = ds.GetGeoTransform()
    left_x = gt[0]
    top_y = gt[3]
    left_x = gt[0] + extent[2]*resolution
    top_y = gt[3] - extent[0]*resolution
    return left_x , top_y

    ## END OF left_top_corner()
###
def latlon_to_rc(ds,lat_lon):

    
    new_cs = osr.SpatialReference()
    new_cs.ImportFromWkt(ds.GetProjectionRef())
    wgs84_wkt = """
    GEOGCS["WGS 84",
        DATUM["WGS_1984",
            SPHEROID["WGS 84",6378137,298.257223563,
                AUTHORITY["EPSG","7030"]],
            AUTHORITY["EPSG","6326"]],
        PRIMEM["Greenwich",0,
            AUTHORITY["EPSG","8901"]],
        UNIT["degree",0.01745329251994328,
            AUTHORITY["EPSG","9122"]],
        AUTHORITY["EPSG","4326"]]"""
    old_cs = osr.SpatialReference()
    old_cs.ImportFromWkt(wgs84_wkt)
    transform = osr.CoordinateTransformation(old_cs,new_cs)
    EN = transform.TransformPoint(lat_lon[1], lat_lon[0]) # swapped 1 and 0 2021-05-15
    #EN = transform.TransformPoint(lat_lon[0], lat_lon[1])
    E = EN[0]
    N = EN[1]
    gt = ds.GetGeoTransform()
    row = int((gt[3]-N)/gt[1])
    col = int((E - gt[0])/gt[1])
    return (row, col)
    ## END OF latlon_to_rc
    
def latlons_to_rcs(ds,lat_lons):

    
    new_cs = osr.SpatialReference()
    new_cs.ImportFromWkt(ds.GetProjectionRef())
    wgs84_wkt = """
    GEOGCS["WGS 84",
        DATUM["WGS_1984",
            SPHEROID["WGS 84",6378137,298.257223563,
                AUTHORITY["EPSG","7030"]],
            AUTHORITY["EPSG","6326"]],
        PRIMEM["Greenwich",0,
            AUTHORITY["EPSG","8901"]],
        UNIT["degree",0.01745329251994328,
            AUTHORITY["EPSG","9122"]],
        AUTHORITY["EPSG","4326"]]"""
    old_cs = osr.SpatialReference()
    old_cs.ImportFromWkt(wgs84_wkt)
    transform = osr.CoordinateTransformation(old_cs,new_cs)
    rows = []
    cols = []
    for lat_lon in lat_lons:
        EN = transform.TransformPoint(lat_lon[1], lat_lon[0]) # swapped 1 and 0 2021-05-15
        #EN = transform.TransformPoint(lat_lon[0], lat_lon[1])
        E = EN[0]
        N = EN[1]
        gt = ds.GetGeoTransform()
        row = int((gt[3]-N)/gt[1])
        col = int((E - gt[0])/gt[1])
        rows.append(row)
        cols.append(col)
    return (rows, cols)
    ## END OF latlon_to_rc

def array_to_raster(ds,array, fname,left_x,top_y,gdal_dtype):
    """Array > Raster
    Save a raster from a C order array.
    :param array: ndarray
    """ 
    dst_filename = fname

    # You need to get those values like you did.
    x_pixels = array.shape[1]  # number of pixels in x
    y_pixels = array.shape[0]  # number of pixels in y
    PIXEL_SIZE = resolution  # size of the pixel...        
    x_min = left_x   
    y_max = top_y  # x_min & y_max are like the "top left" corner.
    wkt_projection = ds.GetProjectionRef()

    driver = gdal.GetDriverByName('GTiff')

    dataset = driver.Create(
        dst_filename,
        x_pixels,
        y_pixels,
        1,
        gdal_dtype, )

    dataset.SetGeoTransform((
        x_min,    # 0
        PIXEL_SIZE,  # 1
        0,                      # 2
        y_max,    # 3
        0,                      # 4
        -PIXEL_SIZE))  

    dataset.SetProjection(wkt_projection)
    dataset.GetRasterBand(1).WriteArray(array)
    dataset.FlushCache()  # Write to disk.
    return dataset, dataset.GetRasterBand(1)  
    ## END OF array_to_raster
    
def easyfix_multiflows(River):

    def generate_juncs_raster(River):

        flow_pointer =  River.flow_pointer
        rows = River.rows
        cols = River.cols
        upwards_pointer = np.full((rows, cols),None, dtype=object)
        juncs_raster = np.zeros((rows,cols))
        for i in np.vstack(np.where(stream_raster > 0)).T:
            try:
                index = D8FlowVals.index(flow_pointer[i[0],i[1]])
                r = i[0] + dY[index]
                c = i[1] + dX[index]
                vals = upwards_pointer[r,c]
                if vals == None:
                    upwards_pointer[r,c] = [inFlowVals[index]]
                else:
                    vals.append(inFlowVals[index])
                    upwards_pointer[r,c] = vals

                juncs_raster[r,c] = len(vals)
            except:
                pass
        return juncs_raster
    
    stream_raster = River.stream_raster
    flow_pointer =  River.flow_pointer
    flow_accum_raster = River.flow_acc
    rows = River.rows
    cols = River.cols
    
    juncs_raster = generate_juncs_raster(River)
    
    def clockwise(point, center):
        shift = (point[0] - center[0], point[1] - center[1])
        n = dYX.index(shift)
        if n == 7:
            return center[0] + dY[0], center[1] + dX[0]
        else:
            return center[0] + dY[n + 1], center[1] + dX[n + 1]

    stream_locations = np.vstack(np.where(stream_raster > 0)).T
    flag = False

    for i in  np.vstack(np.where(juncs_raster > 2)).T:
        row = i[0]
        col = i[1]
        numInFlows = 0
        InFlowLocs = [] 
        for c in range(8):
            x = col + dX[c]
            y = row + dY[c]
            if y < rows and x < cols:
                if (stream_raster[y, x] > 0 and
                        flow_pointer[y, x] == inFlowVals[c]):
                    numInFlows += 1
                    InFlowLocs.append([y, x])
                elif (stream_raster[y, x] > 0 and
                    flow_pointer[row, col] == D8FlowVals[c]):
                    OutFlowCells = [y, x]
        if numInFlows > 2:
            flag = True
            min_area = np.inf
            for i in InFlowLocs:
                
                val = flow_accum_raster[i[0],i[1]]
                if min_area > val:
                    min_area = val
                    mr = i[0]
                    mc =i[1]
            stream_raster[mr,mc] = 0
    return flag #return True if stream_raster has been altered
    ## END OF easyfix_multiflows()




##
class node:
    def __init__(self):
        
        self.ac = None
        self.parent = None
        self.lc = None
        self.rc = None
        self.row = None
        self.col = None
        self.cells = None
        self.category = 'stream'
        self.stream_area_ratio = None
        self.stations = []
        self.void_stations = []
        self.river = None
        
    def add_child(self,data):
        data.parent = self
        if self.rc == None:
            self.rc = data
        else:
            if self.rc.ac > data.ac:
                self.lc = data
            else:
                self.lc = self.rc
                self.rc = data
                
    def NumPost(self):
        if self.category == 'wwtp':
            return int(self.ac * self.stream_area_ratio)
        
        counter = 1
        if self.lc:
            counter += self.lc.NumPost()
        if self.rc:
            counter += self.rc.NumPost()
        return counter



    def pjNum(self):
        counter = 0
        
        if self.parent:
            if self.parent.lc == self:
                counter += self.parent.pjNum()
            else:
                if self.parent.lc:
                    counter += self.parent.pjNum()+self.parent.lc.NumPost()
                else:
                    # leftChild do not exist
                    counter += self.parent.pjNum()
        
        return counter+1



    def pjNum_B(self):
        counter = 0
        
        if self.parent:
            if self.parent.lc == self:
                counter += self.parent.pjNum_B()
            else:
                if self.parent.lc:
                    counter += self.parent.lc.NumPost()+np.max([x.pjNum_B() for x in self.parent.lc.upstreams()])
                else:
                    # leftChild do not exist
                    counter += self.parent.pjNum_B()
        
        return counter+1
    
    def upstreams(self):
        counter = [self]
        if self.lc:
            counter += (self.lc.upstreams())
        if self.rc:
            counter += (self.rc.upstreams())
        return counter
    
    def downstreams(self):
        counter = [self]
        if self.parent:
            counter +=(self.parent.downstreams())
        return counter
    
    def mainstream(self):
        if self.parent:
            if self.parent.lc == self:
                flag = False
            else:
                flag = self.parent.mainstream()
        else:
            flag = True
        return flag
    
    def root(self):
        if self.parent:
            output = self.parent.root()
        else:
            return self
        return output
    
    def calc_ratio(self):
        ratios = []
        for i in self.root().upstreams():
            if i.category == 'stream':
                if len(i.upstreams()) > 5:
                    ratios.append(len(i.upstreams())/float(i.ac))
        return np.mean(ratios)
    
    def pjNum2(self):
        """
        #Normalized Paired Junction order
        #  
        """
        # if not self.river.yscale:
        #     self.river.xscale = (self.river.nodes[0].ac**0.5)
        #     self.river.yscale = self.river.xscale*2.5
        counter = 0.0
        if self.parent:
            if self.parent.lc == self:
                counter = (self.parent.pjNum2())
            else:
                counter = (self.local_outlet().pjNum2()+
                            (self.local_outlet().ac-float(self.ac))/
                            self.root().ac*self.river.yscale)
        
        return counter




        
        return counter
    
    def local_outlet(self): 
        if self.parent:
            if self.parent.lc == self:
                return self
            else:
                return self.parent.local_outlet()
        return self

    def iter_stations(self):
        output = []
        output += self.stations

        if self.lc:
            output += self.lc.iter_stations()
        if self.rc:
            output += self.rc.iter_stations()
        
        return output
    
    def map_stations(self):
        stream_raster = self.river.stream_raster
        import matplotlib.pyplot as plt
        
        plt.matshow(stream_raster)
        xs = []
        ys = []
        for i in self.root().upstreams():
            for stn in i.stations:
                ys.append(stn.snap_loc()[0])
                xs.append(stn.snap_loc()[1])
                #plt.annotate(stn.name,xy=(stn.snap_loc()[1],stn.snap_loc()[0]))
        plt.scatter(xs,ys,color ='yellow')
        plt.show()
    
    def iter_pj(self):
        output = [self]
        if self.lc:
            output += self.lc.iter_pj()
        if self.rc:
            output += self.rc.iter_pj()
        return output
    ## END OF CLASS node
##	
class monitoring_station:

    def __init__(self):
        self.ID = None
        self.loc = None
        self.name = None
        self.stream = None
        self.lat_lon = None
        #self.station_pj =None
        self.source = 'PWQMN'
        self.code= None
        self.lc = None
        self.rc = None
        self.void = False
        self.parent = None
        self.river = None
        
    ##
    def snap_loc(self):
        # point format: [r,c]
        stream_raster = self.river.stream_raster
        flow_pointer = self.river.flow_pointer
        
        point = self.loc
        
        row = point[0]
        col = point[1]

        if row < 0 or col < 0:
            return None
        if row >= self.river.rows or col >= self.river.cols:
            return None

        counter = 0
        while stream_raster[row, col] <= 0:
            if flow_pointer[row, col] == 0:
                return False
            c = D8FlowVals.index(flow_pointer[row, col])
            row = row + dY[c]
            col  = col + dX[c]
            counter += 1
            if counter > 10:
                return False
        return [row, col]
    ##
    def area(self):
        flow_acc= self.river.flow_acc
        r,c = self.snap_loc()
        return flow_acc[r,c]
    ##
    def loc_on_stream(self):
        return self.stream.cells.index(self.snap_loc())
    ## 
    def station_pj(self):
        if not self.void:
            counter = 0
            for river in self.stream.root().iter_pj():
                if river.stations:
                    counter+=len(river.stations)
                if river == self.stream:
                    break
            idx = self.stream.stations[::-1].index(self)
            return counter-idx
        else:
            if self.parent:
                if self.parent.lc == self:
                    return self.parent.station_pj()
                else:
                    return self.lc.station_pj()
            else:
                return 1
    ##
    def lookup(self):
        import webbrowser
        lat,lon = self.lat_lon
        url = 'http://maps.google.com/maps?q=loc:{},{}'.format(lat,lon)
        webbrowser.open_new_tab(url)
    
    def query(self, parm, *masks):
        parm_mask = self.data['PARM'] == parm
        final_mask = True
        for mask in masks:
            final_mask = mask & final_mask
        
        return self.data.loc[parm_mask & final_mask,'RESULT']

    def direct_upstream_stations(self):
        stream = self.stream
        waitlist = [stream]
        output = []

        while waitlist:
            cursor = waitlist[0]
            waitlist.remove(cursor)

            if cursor.lc:
                if cursor.lc.stations:
                    output.append(cursor.lc.stations)
                else:
                    waitlist.append(cursor.lc)
            
            elif cursor.rc:
                if cursor.rc.stations:
                    output.append(cursor.lc.stations)
                else:
                    waitlist.append(cursor.rc)
        return(output)
    ## END OF CLASS monitoring_station
##

class river:

    def __init__(self):
        self.name = None
        self.flow_pointer = None
        self.stream_raster = None
        self.flow_acc = None
        self.rows = None
        self.cols = None
        self.lat_lon = None
        self.database = None
        self.nodes = None
        self.outlet = None
        self.stations = None
        self.extent = None
        self.alerts = {}
        self.failed_stations = []
        self.yscale = None
        self.xscale = None
   

    def generate_nodes(self):
        """
        This function does 3 things: 1. Read the 

        Args:
            self (pj.river object)

        Changes:
            - 

        Returns:
            None

        Requirements:

            - 

        UsedIn:
            - 
        """
        #establish a paired junction network

        stream_raster = self.stream_raster
        flow_pointer = self.flow_pointer
        flow_accum_raster = self.flow_acc
        rows = self.rows
        cols = self.cols

        rivermouth = self.outlet
        root = node()
        root.river = self
        root.ac = flow_accum_raster[rivermouth[0],rivermouth[1]]
        root.row, root.col = rivermouth
        waitlist = [root]

        ID = 0
        counter = 0
        
        for i in waitlist:
            counter +=1
            row = i.row
            col = i.col
            inlets = [i]
            numIn =1
            cellsList = [[row,col]]
            while numIn == 1:      
                inlets = []
                
                numIn = 0
                for c in range(8):
                    y = row + dY[c]
                    x = col + dX[c]
                    if y< rows and x < cols:
                        if (stream_raster[y,x]>0 and 
                            flow_pointer[y,x] == inFlowVals[c]):
                            
                            inlets.append([y,x])
                            numIn +=1
                if numIn == 1:
                    row = inlets[0][0]
                    col = inlets[0][1]
                    cellsList.append([row,col])
                    
            i.cells = cellsList

            for inlet in inlets:
                ID+=1
                tmpNode = node()
                tmpNode.row = inlet[0]
                tmpNode.col = inlet[1]
                tmpNode.ac = flow_accum_raster[inlet[0],inlet[1]]
                tmpNode.river = self
                if stream_raster[inlet[0],inlet[1]] == 2:
                    tmpNode.category = 'wwtp'
                else:
                    tmpNode.category = 'stream'
                waitlist.append(tmpNode)
                i.add_child(tmpNode)
        self.nodes = waitlist

    def attach_stations(self):


        """
        This function does 3 things: 1. Read the 'PWQMN_STATIONS.csv' and create 
        monitoring_station objects 2. Spatially join station objects to nodes objects.
        3. If there are mutiple stations on one node, re-order them in sequence.

        Args:
            self (pj.river object)

        Changes:
            - 

        Returns:
            None

        Requirements:

            - 

        UsedIn:
            - 
        """

        stream_raster = self.stream_raster
        flow_pointer = self.flow_pointer
        flow_accum_raster = self.flow_acc
        rows = self.rows
        cols = self.cols
        alert_message = pd.DataFrame()

        PWQMN = pd.read_csv('PWQMN_STATIONS.csv', dtype = object)
        PWQMN['LONGITUDE'] = pd.to_numeric(PWQMN['LONGITUDE'])
        PWQMN['LATITUDE'] = pd.to_numeric(PWQMN['LATITUDE'])
        
        lat_lons = np.array([PWQMN['LONGITUDE'],PWQMN['LATITUDE']]).T
        rs, cs = latlons_to_rcs(self.ds,lat_lons)
        PWQMN['row'] = rs
        PWQMN['col'] = cs
        
        output = []
        for index, row in PWQMN.iterrows():
            a = [row['row'],row['col']]
            extent = self.extent
            loc = [a[0]-extent[0], a[1]-extent[2]]
            #pdb.set_trace()
            if loc[0] < rows and loc[0] > 0:
                if loc[1] <cols and loc[1] > 0:
                    if flow_pointer[loc[0],loc[1]] > 0:
                        out = monitoring_station()
                        out.loc = loc
                        out.ID = row['STATION']
                        out.river = self
                        if out.ID.startswith('040013_LD_'):
                            out.source = 'LONDON'
                        out.name = [row['NAME'],row["LOCATION"]]
                        out.lat_lon = [row["LATITUDE"],row["LONGITUDE"]]
                        output.append(out)
        
        
        ## assign 'stream' attribute to each station object
        station_code = pd.read_csv('database/Crawley_Stations_Code.csv'
                        ,index_col = 'ID')['Code'].to_dict()
        for stn in output:
            if stn.ID in list(station_code.keys()):
                stn.code = station_code[stn.ID]

            for river in self.nodes:
                if stn.snap_loc() in river.cells:
                    stn.stream = river
                    river.stations += [stn]
            if stn.stream == None:
                self.failed_stations.append(stn)
                arow = pd.DataFrame([[stn.ID,stn.code,stn.name[0],stn.name[1]]],
                    columns=['STATION','CODE','NAME','DESCRIPTION'])
                #alert_message = alert_message.append(arow,ignore_index=True) #temporary block

        #self.alerts['failed_stations']=alert_message # temporary block

        ## re-order stations on a same stream segment; This 
        for river in self.nodes[0].root().iter_pj():
            if len(river.stations) > 1:
                idx_dict = {}
                for stn in river.stations:
                    idx= river.cells.index(stn.snap_loc())*100
                    while idx in idx_dict:
                        idx += 1       
                    idx_dict[idx] = stn
                list(idx_dict.keys()).sort()
                river.stations = []
                for idx in sorted(idx_dict.keys()):
                    river.stations.append(idx_dict[idx])
        self.stations = self.nodes[0].root().iter_stations()
    
    def merge_duplicates(self):
        
        """
        Identify stations too close to each other (5 cells/150m). Merge the
        stations by unifying their ID in database. Station objects merged 
        will be removed from river.stations.

        Note: remeber to run split_database() after using this function

        Args:
            self (pj.river object)

        Changes:
            - Station objects merged will be removed from the river.stations
            list
            - In database, records for deleted stations with be modified to 
                the new station's ID

        Returns:
            None

        Requirements:

            - river.nodes (lists): with nodes objects, which has node.stations 
            - river.database

        UsedIn:
            - river.quick function
        """
        database = self.database 
        alert_message = pd.DataFrame(columns=['STATIONS'])

        for nd in self.nodes:
            if len(nd.stations)>1:
                
                for stn in nd.stations:
                    removed_stations = []
                    flag = False
                    for stn2 in list(set(nd.stations)-set([stn])):
                        if stn !=  stn2:
                            if abs(nd.cells.index(stn.snap_loc())
                                - nd.cells.index(stn2.snap_loc())) <10:
                                if stn.source == stn2.source:
                                    mask = database['STATION'] == stn2.ID
                                    database.loc[mask,'STATION']= stn.ID
                                    nd.stations.remove(stn2)
                                    self.stations.remove(stn2)
                                    if stn2.code and (not stn.code):
                                        stn.code = stn2.code
                                    removed_stations.append(stn2)
                                    flag = True
                                else:
                                    pass
                            else:
                                # add pseduo stream
                                pass 
                        else:
                            pass
                    if flag:
                        alert_message.loc[stn.ID,'STATIONS'] = ','.join([st.ID for st in removed_stations])
            else:
                pass
                                    
        self.alerts['merges']= alert_message 
    
    
    def stations_dict(self):
        '''
        This function will build a dictionary for all station objects
        in self.stations. 
        
        Args:
            self (pj.river object)
        Changes:
            None
        Returns:
            stations_dict (dict): a dictionary contains all stations object
                where station ID are keys
        '''
        tmp_dict = {}
        for stn in self.stations:
            tmp_dict[stn.ID] = stn
        return tmp_dict
        # END of stations_dict()


    def load_raster(self, threshold=5):
        '''
        This function read rasters (flow_pointer,flow_acc) for this river object. 
        
        Args:
            self (pj.river object)
        Changes:
            self.database (pd.DataFrame)

        '''

        ds = gdal.Open('raster_folder/EnhancedFlowDirection_ON2-edit.tif')
        self.ds = ds
        cols = ds.RasterXSize
        rows = ds.RasterYSize
        flow_pointer_ON = ds.ReadAsArray(0, 0, cols, rows).astype(int)
        RC = latlon_to_rc(self.ds, self.lat_lon)
        self.RC_ON = RC
        ## external function latlon_to_rc()

        self.snap_outlet()
        RC = self.RC_ON
        history = pd.read_csv('cache/history.csv')

        if RC not in np.array(history):

            ## trace catchment
            tmpGrid = np.full((rows, cols), -99, dtype=int)
            list_y = []
            list_x = []
            waitlist = [RC]
            for cell in waitlist:
                row = cell[0]
                col = cell[1]
                val = flow_pointer_ON[row, col]
                tmpGrid[row, col] = val
                for c in range(8):
                    try:
                        x = col + dX[c]
                        y = row + dY[c]
                        if flow_pointer_ON[y, x] == inFlowVals[c]:
                            waitlist.append([y, x])
                            list_y.append(y)
                            list_x.append(x)
                    except:
                        pass
            max_r = max(list_y)
            min_r = min(list_y)
            max_c = max(list_x)
            min_c = min(list_x)
            extent = [min_r, max_r, min_c, max_c]
            flow_pointer = tmpGrid[min_r:max_r+1, min_c:max_c+1]
            flow_acc = Flow_Acc(flow_pointer)

            ## external function Flow_Acc()
            left,top = left_top_corner(self.ds, extent)
            array_to_raster(self.ds, flow_pointer,'cache/pointer/{}_{}.tif'.format(RC[0],RC[1]),
                            left, top, gdal.GDT_Int16)
            array_to_raster(self.ds, flow_acc,'cache/FlowAcc/{}_{}.tif'.format(RC[0],RC[1]),
                            left, top, gdal.GDT_Float32)
            pd.to_pickle(extent,'cache/extent/extent_{}_{}.pickle'.format(RC[0],RC[1]))

            df_in = pd.DataFrame(data={'row': RC[0], 'col': RC[1]}, index=[0])[['row','col']]
            history.append(df_in).to_csv('cache/history.csv', index=False)

        else:
            flow_pointer = read_raster('cache/pointer/{}_{}.tif'.format(RC[0],RC[1]), int)
            flow_acc = read_raster('cache/FlowAcc/{}_{}.tif'.format(RC[0],RC[1]), float)
            extent = pd.read_pickle('cache/extent/extent_{}_{}.pickle'.format(RC[0],RC[1]))
            
        
        
        self.extent = extent
        self.flow_pointer = flow_pointer
        self.outlet = [RC[0]-extent[0], RC[1]-extent[2]]
        self.rows = self.flow_pointer.shape[0]
        self.cols = self.flow_pointer.shape[1]
        self.flow_acc = flow_acc

        stream_raster = np.array(self.flow_acc)
        stream_raster[stream_raster < threshold ] = -1
        stream_raster[stream_raster >= threshold ] = 1
        stream_raster = stream_raster.astype(int)

        self.stream_raster = stream_raster
        
        ## fix mulit inflows issue
        easyfix_multiflows(self)

        ## END OF load_raster()

    def load_database(self):
        '''
        This function select relevant records from PWQMN database ('PWQMN_SQLite2.db'). 
        Note self.stations much be valid in order to select all records.
        
        Args:
            self (pj.river object)
        Changes:
            self.database (pd.DataFrame)
        Returns:
            None
        
        '''
        import sqlite3
        station_ID = tuple([stn.ID for stn in self.stations])
        conn = sqlite3.connect('PWQMN_SQLite2.db')
        sql = 'SELECT * FROM PWQMN_database WHERE STATION IN {}'.format(station_ID)
        df = pd.read_sql(sql, conn,parse_dates='DATE')
        #df['TIME'] = pd.to_timedelta(df['TIME'])
        self.database = df

    def split_database(self):
        """
        This function split self.database accroding to station ID and
        assign each sub-database to each station (as station.data).
        Note: If there is a station with no data, it will be assigned a
        DataFrame with all columns but not no row (empty header).

        Args:
            self (pj.river object)
        Changes:
            station.data (pd.DataFrame): It will assign sub-database to
                each station object
        Returns:
            None
        """
        gp_database = self.database.groupby('STATION')
        for stn in self.stations:
            if stn.ID in list(gp_database.indices.keys()):
                stn.data = gp_database.get_group(stn.ID)
            else:
                stn.data = self.database.head(0).copy()

            

    def back_data_up(self):

        self.database_backup = self.database.copy()
        for stn in self.stations:
            stn.data_backup = stn.data.copy()

    def recover_data(self):
        self.database = self.database_backup.copy()
        for stn in self.stations:
            stn.data = stn.data_backup.copy()
        ## END OF recover_data

    def add_extra_fileds(self):

        database = self.database
        stn_pj_dict = dict([[stn.ID, stn.station_pj()] for stn in self.stations])
        database['PJ_ORDER'] = database['STATION'].map(stn_pj_dict)
        database['STN_DT'] = database['STATION'] + database['DATE'].dt.strftime(' %Y-%m-%d')


    def xyscale(self):
        if not self.yscale:
            self.xscale = (self.nodes[0].ac**0.5)
            self.yscale = self.xscale*2.5

    def quick(self,threshold=8.0):
        """
        harness function
        threshold: channalizing threshold (KM^2)
        """

        print('loading rasters.....')
        self.load_raster(threshold)
        print('generate nodes.....')
        self.generate_nodes()
        self.attach_stations()
        print('loading database.....')
        self.load_database()
        print('merge ducpliates.....')
        self.merge_duplicates()
        self.split_database()
        self.add_extra_fileds()
        self.back_data_up()
        self.xyscale()
        print('completed')
        #END of quick()
        
    def snap_outlet(self):
        
        """
        This function are used to snap the ............

        Args:
            self (pj.river object)
        Changes:


        Returns:
            None
        Used in:
            load_raster()
        """
        
        rc = self.RC_ON 
        StreamGRID = read_raster('raster_folder/StreamGRID.tif',bool)
        flow_pointer = read_raster('raster_folder/EnhancedFlowDirection_ON2-edit.tif',int)
        counter = 0
        while not StreamGRID[rc[0], rc[1]] and counter <10:
            counter +=1
            index = D8FlowVals.index(flow_pointer[rc[0], rc[1]])
            row = rc[0] + dY[index]
            col = rc[1] + dX[index]
            rc = [row, col]
        self.RC_ON = rc

        if not StreamGRID[rc[0],rc[1]]:
            print('failed to snap rivermouth in 10 cells distance')
            raise

    def query(self, parm, *masks):
        parm_mask = self.database['PARM'] == parm
        final_mask = True
        for mask in masks:
            final_mask = mask & final_mask
        
        return self.database.loc[parm_mask & final_mask,'RESULT']
        ## END OF query


    def create_pjRaster(self):

        flow_pointer = self.flow_pointer
        rows = self.rows
        cols = self.cols
        flow_acc = self.flow_acc
        mx = np.count_nonzero(flow_pointer>0)

        pjGrid = np.zeros((rows, cols),dtype=int)
        waitlist = [self.outlet]
        counter = 0

        while len(waitlist)>0: 

            cursor = waitlist[0]
            #idx = D8FlowVals.index(flow_pointer[cursor[0], cursor[1]])
            xs = []
            ys = []
            pjGrid[cursor[0], cursor[1]] = counter
            counter+=1
            waitlist = waitlist[1:]

            for c in range(8):
                row = cursor[0] + dY[c]
                col = cursor[1] + dX[c]
                if row < rows and col < cols:
                    if flow_pointer[row, col] == inFlowVals[c]:
                        xs.append([row, col])
                        ys.append(flow_acc[row,col])

            arr = np.array(list(zip(xs,ys)))
            if arr.size == 0:
                continue
            arr = arr[arr[:, 1].argsort()]
            waitlist = arr[:,0].tolist() +waitlist

            if counter > (mx+100):
                print('interrupted')
                break
            if counter%500000 == 0:
                print((round(float(counter)/(mx+10.)*100,1),'%'))

        self.pj_raster = pjGrid
        print('completed')


    def __getitem__(self, args):
        if type(args) is not tuple:
            args =[args]

        final_mask = True
        for arg in args:
            mask = True
            if arg in self.database['PARM'].unique():
                mask = self.database['PARM'] == arg
            elif arg in self.database['STATION'].unique():
                mask = self.database['STATION'] == arg
            elif '-' in arg:
                mask = self.database['DATE'] == pd.to_datetime(arg)
            final_mask = mask & final_mask
        return self.database[final_mask]

    def pivot_table(self,*args):
        """
        return a pivot table based on PARM
        """
        mask = self.database['PARM'].isin(args)
        df = self.database[mask]
        pivot = df.pivot_table(index=['STN_DT'], 
                                columns='PARM', 
                                values='RESULT',
                                dropna=True)
        pivot = pivot.dropna(axis=0, how='any')
        pivot['STATION'] = pivot.index.map(lambda x: x.split(' ')[0])
        pivot['DATE'] = pivot.index.map(lambda x: pd.to_datetime(x.split(' ')[1]))
        return pivot

## End of class river

class interpolator:

    def __init__(self,River):
        self.river = River
        River.interpolator = self
        self.parm = None
        self.date1 = None
        self.date2 = None

        self.stat_func = None
        self.colormap = None
        self.gamma = None
        self.min_sample = None
        self.zmin= None
        self.zmax = None
        self.disable_kml_export = False
        self.MonthFilter = [1,2,3,4,5,6,7,8,9,10,11,12]

    def SetValues(self, parm='COND25', year1 =1980, year2 =1995,
                stat_func=np.median, colormap=matplotlib.cm.rainbow,
                gamma=1, min_sample=20):
        """
        Set the major parameters before kriging
        
        """
        self.parm = parm
        self.date1 = pd.to_datetime('{}-1-1'.format(year1))
        self.date2 = pd.to_datetime('{}-12-31'.format(year2))
        self.stat_func = stat_func
        self.colormap = colormap
        self.gamma = gamma
        self.min_sample = min_sample

    def SetPercentile(self,perct):
        """
        """
        self.stat_func = lambda x: np.percentile(x, perct)
        self.percentile = perct
    
    def kriging(self):
        
        from pykrige.ok import OrdinaryKriging
        import pykrige.kriging_tools as kt
        import copy
        self.stations_copy = []
        self.xscale = self.river.xscale
        self.yscale = self.river.yscale

        dt1 = self.date1
        dt2 = self.date2

        print(('minimum {} samples'.format(self.min_sample)))

        from collections import Counter
        database = self.river.database
        b = Counter(database[database['PARM'] == self.parm]['UNITS'])
        self.unit = b.most_common(1)[0][0]

        xs,ys,zs,sampled_stations = [],[],[],[]
        sampled_data = []
        cts = []
        for stn in self.river.stations:
            mask1 = stn.data['PARM'] == self.parm
            mask2 = stn.data['DATE'] > dt1
            mask3 = stn.data['DATE'] < dt2
            mask4 = stn.data['DATE'].dt.month.isin(self.MonthFilter)
            result = stn.data[mask1 & mask2 & mask3 & mask4]['RESULT']

            stn.sampled = False
            if result.size >= self.min_sample:
                sampled_stations.append(stn)
                sampled_data.append(result)
                xs.append(np.sqrt(stn.area()))
                ys.append(stn.stream.pjNum2())
                zs.append(self.stat_func(result))
                cts.append(result.size)
                #ct.append(result.size)
                stn.sampled = True
                stn.tmp_zs = [round(n,3) for n in result.tolist()]
                stn.tmp_zs_dt = stn.data.loc[result.index]['DATE']
                stn.tmp_z = self.stat_func(result)
                self.stations_copy.append(copy.copy(stn))

        self.max_x = np.sqrt(self.river.nodes[0].ac)
        self.zmin = np.min(zs)
        self.zmax = np.max(zs)
        self.sampled_stations = sampled_stations
        self.sampled_data = sampled_data
        self.xs = xs
        self.ys = ys
        self.cts = cts

        xs = [x/self.max_x*self.xscale for x in xs] #Normalize X

        data = np.array(list(zip(xs,ys,zs)))
        data = data[~np.isnan(data).any(axis=1)]

        self.zs = data[:,2]
        self.zdata = data
        
        if data.size == 0:
            print('no data available')

        gridx = np.arange(0.,self.xscale*1.1, 1)
        gridy = np.arange(0., self.yscale*1.1, 1)

        OK = OrdinaryKriging(data[:, 0], data[:, 1], data[:, 2], 
                         variogram_model='exponential',
                         variogram_parameters=[30,50,2],
                         verbose=True, enable_plotting=False)
        
        z, ss = OK.execute('grid', gridx, gridy)
        self.OK = OK
        self.z_matrix = z



    def inverse_trans(self, display=False, save=False,format='kml'):
        import copy

        display_mat = np.empty((self.river.rows, self.river.cols))
        display_mat[:] = np.NAN
        self.nodes_copy = [] # create an unique node objects copy for interpolator

        for stream in self.river.nodes:
            self.current_node = stream
            x = int(np.sqrt(stream.ac/self.max_x*self.xscale) + 0.5)
            y = int(stream.pjNum2() + 0.5)
            setattr(stream, self.parm, self.z_matrix[y,x])
            stream.tmp_z = self.z_matrix[y,x]
            self.nodes_copy.append(copy.copy(stream)) #
            
            for cell in stream.cells:
                display_mat[cell[0],cell[1]] = getattr(stream, self.parm)
        
        if save:
            if format == 'shp':
                self.write_shapefile()
            elif format == 'kml':
                self.write_kml()

            elif format == 'time_lapse':
                self.write_time_lapse()
        if display:
            plt.matshow(display_mat)
            plt.show()
                
    def network_interp(self, save=False,format='kml'):
        import copy
        self.nodes_copy = [] # create an unique node objects copy for interpolator

        #set pv for each station
        for nd in self.river.nodes:
            nd.pv = None
            nd.av = None
        for stn in self.river.stations:
            if stn.sampled:
                if not stn.stream.pv:
                    stn.stream.pv = stn.tmp_z
                elif stn.source == "LONDON":
                    stn.stream.pv = stn.tmp_z

        #set av for each station
        for nd in self.river.nodes:
            counter = 0
            if nd.pv:
                load = nd.pv * nd.ac
                ac = nd.ac
                waitlist = [nd]
                while waitlist:
                    cd = waitlist[0]
                    waitlist.remove(cd)
                    if cd.lc:
                        if cd.lc.pv:
                            load = load - cd.lc.pv * cd.lc.ac
                            ac = ac-cd.lc.ac
                        else:
                            waitlist.append(cd.lc)
                            
                    if cd.rc:
                        if cd.rc.pv:
                            load = load - cd.rc.pv * cd.rc.ac
                            ac = ac-cd.rc.ac
                            
                        else:
                            waitlist.append(cd.rc)
                    counter +=1
                nd.av = load/ac

        #set pv, av for headwaters, av for intermediate streams
        for nd in self.river.nodes: 
            if not nd.pv:
                if nd.parent:       
                    tmp_nd = nd
                    flag = True
                    while not tmp_nd.parent.pv:
                        tmp_nd = tmp_nd.parent
                        if not tmp_nd.parent:
                            #print 'nan'
                            flag = False
                            break
                    if flag:
                        if not(nd.lc or nd.rc):
                            nd.pv = tmp_nd.parent.av
                            nd.av = nd.pv
                        else:
                            nd.av = tmp_nd.parent.av

        #interation:
        counter = 1
        counter0 =0
        while counter>counter0:
            counter0 = counter
            for nd in self.river.nodes:
                if not nd.pv:
                    if nd.lc and nd.rc:
                        if nd.lc.pv and nd.rc.pv and nd.av:
                            counter+=1
                            nd.pv = (nd.lc.pv* nd.lc.ac
                                    +nd.rc.pv* nd.rc.ac
                                    +(nd.ac-nd.lc.ac-nd.rc.ac)*nd.av )/nd.ac


        for stream in self.river.nodes:
            if stream.pv:
                setattr(stream, self.parm, stream.pv)
                stream.tmp_z = stream.pv
                cp = copy.copy(stream)
                cp.original = stream
                stream.cp = cp
                self.nodes_copy.append(cp)


        self.zmin = np.percentile([nd.pv for nd in self.nodes_copy],10)
        self.zmin = 0 #ad hoc
        self.zmax = np.percentile([nd.pv for nd in self.nodes_copy],99)

        if save:
            if format == 'shp':
                self.write_shapefile()
            elif format == 'kml':
                self.write_kml()

            elif format == 'time_lapse':
                self.write_time_lapse()
        #Clean up 

        for nd in self.river.nodes:
            nd.pv = None
            nd.av = None

    def write_shapefile(self):
        
        import fiona
        import os
        from shapely.geometry import Point, mapping, LineString
        import uuid

        ## Create a directory

        ## Create a directory
        if not os.path.exists("shp/{}".format(self.parm)):
            os.makedirs("shp/{}".format(self.parm))


        fname = '{}_{}_{}-{}'.format(self.river.name.split()[0],
            self.parm, str(self.date1.year)[-2:], str(self.date2.year)[-2:])
            
        if self.percentile:
            fname = "{}pct_".format(self.percentile) + fname

        unique_folder = uuid.uuid4().hex[:6].upper()
        dirr = "shp/{}/{}_{}".format(self.parm, unique_folder, fname)
        os.makedirs(dirr)


        ## Write Streams Shapefile (line)
        schema = { 'geometry': 'LineString', 'properties': { self.parm: 'float' } }
        left, top = left_top_corner(self.river.ds, self.river.extent)
        prj = self.river.ds.GetProjection()



        with fiona.collection("{}/STREAM_{}.shp".format(dirr, fname), 
            "w", "ESRI Shapefile", schema, crs=prj) as output1:

            for stream in self.nodes_copy:
                if len(stream.cells)>1:
                    points = [[left +cell[1]*30+15, top-cell[0]*30-15] for cell in stream.cells]
                    points = points[0::3]+[points[-1]] # simplify lines
                    line = LineString(points)
                elif len(stream.cells) == 1:
                    points = [[left+stream.cells[0][1]*30, top - stream.cells[0][0]*30]]*2
                    line = LineString(points)
                else:
                    print('an error occured when creating stream polyline, skipped')
                    continue
        
                output1.write({
                    'properties': {
                        self.parm: getattr(stream, self.parm)
                    },
                    'geometry': mapping(line)
                })

        ## Write Stations Shapefile (point)
        schema = { 'geometry': 'Point', 'properties': { self.parm: 'float' } }
        wkt= ('GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84"'
        ',6378137,298.257223563,AUTHORITY["EPSG","7030"]],'
        'AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,'
        'AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,'
        'AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]]')


        with fiona.collection("{}/STN_{}.shp".format(dirr, fname),
            "w", "ESRI Shapefile", schema, crs=wkt) as output2:
            
            for stn in self.river.stations:
                if stn.sampled == True:
                    point = Point([stn.lat_lon[1], stn.lat_lon[0]])
                    output2.write({
                        'properties': {
                            self.parm: getattr(stn, 'tmp_z'),
                        },
                        'geometry': mapping(point)
                    })


    def write_kml(self):

        from pykml.factory import KML_ElementMaker as KML
        from lxml import etree
        import osr
        import os
        import uuid

        ## Create a directory
        if not os.path.exists("kml/{}".format(self.parm)):
            os.makedirs("kml/{}".format(self.parm))


        fname = '{}_{}_{}-{}'.format(self.river.name.split()[0],
            self.parm, str(self.date1.year)[-2:], str(self.date2.year)[-2:])

        if self.percentile:
            fname = "{}pct_".format(self.percentile) + fname

        unique_folder = uuid.uuid4().hex[:6].upper()
        dirr = "kml/{}/{}_{}".format(self.parm, unique_folder, fname)
        os.makedirs(dirr)
        os.makedirs("{}/plots".format(dirr))

        ## Construct KML
        name_object = KML.name("stream")
        Style1 = KML.Style(KML.IconStyle(KML.scale(0.7),KML.Icon(KML.href(
                        'http://maps.google.com/mapfiles/kml/pal3/icon61.png'))),
                        KML.LabelStyle(KML.scale(0.8)),id = "s_ylw-pushpin")

        Style2 = KML.Style(KML.IconStyle(KML.scale(0.7),KML.Icon(KML.href(
                        'http://maps.google.com/mapfiles/kml/pal3/icon61.png'
                        ))),id = "s_ylw-pushpin_hl")
        StyleMap = KML.StyleMap(
            KML.Pair(KML.key('normal'),KML.styleUrl('#s_ylw-pushpin')),
            KML.Pair(KML.key('highlight'),KML.styleUrl('#s_ylw-pushpin')),
            id = "m_ylw-pushpin")
        visb = KML.visibility('0')
        name_tag = KML.name('Stations')
        TimeSpan = KML.TimeSpan(KML.begin('{:%Y-%m-%d}'.format(self.date1)),
                        KML.end('{:%Y-%m-%d}'.format(self.date2)))
        self.TimeSpan = TimeSpan

        ## POINTS
        PM = {}
        folder1 = KML.Folder(name_tag, Style1,Style2,StyleMap) #visb
        StrTimeWindow = self.date1.strftime('%Y/%m/%d') + self.date2.strftime('-%Y/%m/%d')
        
        for stn in self.stations_copy:
            if stn.sampled:
                # plot data
                plt.plot(stn.tmp_zs_dt, stn.tmp_zs,marker='.',linestyle='None')
                imgURL = 'plots/{}.png'.format(uuid.uuid4().hex[:5].upper())
                imgtag = '<img src="{}" />'.format(imgURL)

                plt.xlim(self.date1,self.date2)
                plt.gcf().set_size_inches(8,6)
                plt.gca().margins(tight=True)
                plt.savefig(dirr+"/"+imgURL, bbox_inches='tight')
                plt.clf()


                yx = str(stn.lat_lon[1]) +',' +str(stn.lat_lon[0])
                vz = float("{0:.3g}".format(stn.tmp_z)) #visual display value (rounded)
                if vz>=10:
                    vz = round(vz,1)
                data_tag = KML.ExtendedData(
                        KML.Data(KML.value(stn.ID), name='Station ID'),
                        KML.Data(KML.value(stn.name[0]), name='Stream'),
                        KML.Data(KML.value(stn.name[1]), name='Site'),
                        KML.Data(KML.value(StrTimeWindow), name='Time Window'),
                        KML.Data(KML.value(str(len(stn.tmp_zs))), name='Sample Size'),
                        KML.Data(KML.value(imgtag), name='Plot')
                                            )

                PM[stn.ID] = KML.Placemark(KML.name(vz),data_tag,
                                #KML.visibility('0'),
                                KML.Point(KML.coordinates(yx)),
                                KML.styleUrl('#m_ylw-pushpin'))
                folder1.append(PM[stn.ID])

        ## COLORMAP

        self.colormap.set_gamma(self.gamma)

        def convert_kml_color(val):
            """
            Convert parm value to an color on self.colormap
            """
            z = (val-self.zmin)/(self.zmax-self.zmin)
            rgb = tuple((np.array(self.colormap(z))*255).round(3)[:3])
            hex6 = '%02x%02x%02x' % rgb
            return 'FF' + hex6[4:6]+hex6[2:4] +hex6[0:2]

        ## LINE
        def rc_to_latlon(r_c_list):
            """
            Convert a list of [row, col] to [lat, lon]
            """
            ds = self.river.ds
            gt = ds.GetGeoTransform()
            output = []
        
            old_cs = osr.SpatialReference()
            old_cs.ImportFromWkt(ds.GetProjectionRef())
            wgs84_wkt = """
            GEOGCS["WGS 84",
                DATUM["WGS_1984",
                    SPHEROID["WGS 84",6378137,298.257223563,
                        AUTHORITY["EPSG","7030"]],
                    AUTHORITY["EPSG","6326"]],
                PRIMEM["Greenwich",0,
                    AUTHORITY["EPSG","8901"]],
                UNIT["degree",0.01745329251994328,
                    AUTHORITY["EPSG","9122"]],
                AUTHORITY["EPSG","4326"]]"""
            new_cs = osr.SpatialReference()
            new_cs.ImportFromWkt(wgs84_wkt)
            transform = osr.CoordinateTransformation(old_cs,new_cs)
            
            for r_c in r_c_list:
                row =self.river.extent[0] + r_c[0]
                col =self.river.extent[2] + r_c[1]
                y = gt[3] - row*30
                x = gt[0] + col*30
                LATLON = transform.TransformPoint(x, y)[::-1][1:]
                LATLON = (round(LATLON[0],4), round(LATLON[1],4)) # round to reduce kml size
                output.append(LATLON)
                
            return output
            ## END OF rc_to_latlon

        for nd in self.nodes_copy:
            if nd.parent:
                tail = [nd.parent.cells[-1]] # this close the joints.
            else:
                tail = []
            nd.cells_latlon = rc_to_latlon(tail + nd.cells)

        ## Construct line KML
        name_tag = KML.name('Streams')
        
        folder2 = KML.Folder(name_tag)
        PM2 ={}
        counter = 0

        for nd in self.nodes_copy:
            counter += 1
             
            xy_string= ''
            for cell in nd.cells_latlon[0::4] + [nd.cells_latlon[-1]]:
                txt = str(cell[1]) +','+ str(cell[0])+' '
                xy_string += txt

            PM2[counter] = KML.Placemark(KML.name("stream"),
                                KML.ExtendedData(KML.Data(KML.value(round(nd.tmp_z,4)), name=self.parm),
                                                 KML.Data(KML.value(str(round(nd.ac,3))+' KM2'),name='Catchment')),
                                KML.Style(KML.LineStyle(
                                    KML.color(convert_kml_color(nd.tmp_z)),
                                                    KML.width('4.58')),
                                        KML.PolyStyle(KML.fill('0'))),
                                KML.LineString(KML.altitudeMode('clampToGround'),
                                            KML.coordinates(xy_string)))
            folder2.append(PM2[counter])


        

        ## Construct Legend Overlay
        ScrOv = KML.ScreenOverlay(KML.name('legend'),
                                KML.Icon(KML.href(fname + ".png")),
                                KML.overlayXY(x="0.005",y="0.005",
                                            xunits="fraction",yunits="fraction"),
                                KML.screenXY(x="0.005",y="0.005",
                                            xunits="fraction",yunits="fraction"),
                                KML.rotationXY(x="0",y="0",xunits="fraction",
                                            yunits="fraction"),
                                KML.size(x="0",y="0",xunits="fraction",
                                            yunits="fraction"),
                                )
        

        ## Assign kml strings to self object
        self.kml_folder1 = folder1
        self.kml_folder2 = folder2
        self.kml_ScrOv = ScrOv

        if self.disable_kml_export:
            self.disable_kml_export = False
            return

        ## Plot A Legend.png

        

        fig = plt.figure(figsize=(6, 1))
        ax1 = fig.add_axes([0.05, 0.80, 0.9, 0.15])
        ax1.tick_params(axis='x', colors='white')
        norm = matplotlib.colors.Normalize(self.zmin, self.zmax)
        cb1 = matplotlib.colorbar.ColorbarBase(ax1, cmap=self.colormap,
                                norm=norm,
                                orientation='horizontal')
        cb1.set_label(self.parm + ' / ' + self.unit,color='white')
        plt.savefig("{}/{}.png".format(dirr, fname), transparent=True)
        plt.clf()
        plt.close()
        
        ## Construct output KML
        kml_tag = KML.kml(KML.Document(ScrOv,folder1,folder2))
        doc = etree.tostring(kml_tag,pretty_print=True)

        kml_file = open("{}/{}.kml".format(dirr, fname), "w")
        kml_file.write(doc)
        kml_file.close()

## End of class interpolator()


class Province(river):

    def __init__(self):
        river.__init__(self)
        self.ds = gdal.Open('raster_folder/EnhancedFlowDirection_ON2-edit.tif')
        self.ds2 = gdal.Open('raster_folder/FlowACC_ON2_float2.tif')


    def pv_load_raster(self):

        self.flow_pointer = self.ds.ReadAsArray(0, 0, self.ds.RasterXSize,
                            self.ds.RasterYSize).astype(int)
        self.flow_accum_raster = self.ds2.ReadAsArray(0, 0, self.ds2.RasterXSize,
                                 self.ds2.RasterYSize).astype(float)

        self.stream_raster = (self.flow_accum_raster>5).astype(int)

        self.extent = [0, self.ds.RasterYSize-1, 0, self.ds.RasterXSize-1]
        self.rows = self.ds.RasterYSize
        self.cols = self.ds.RasterXSize

    def pv_load_stations(self):
                #attach Monitoring Stations
        PWQMN = pd.read_csv('PWQMN_STATIONS.csv', dtype = object)
        PWQMN['LONGITUDE'] = pd.to_numeric(PWQMN['LONGITUDE'])
        PWQMN['LATITUDE'] = pd.to_numeric(PWQMN['LATITUDE'])
        
        output = []

        for index, row in PWQMN.iterrows():
            a = latlon_to_rc(self.ds, [row["LONGITUDE"], row["LATITUDE"]])
            extent = self.extent
            loc = [a[0]-extent[0], a[1]-extent[2]]

            out = monitoring_station()
            out.loc = loc
            out.ID = row['STATION']
            out.river = self
            if out.ID.startswith('040013_LD_'):
                out.source = 'LONDON'
            out.name = [row['NAME'],row["LOCATION"]]
            out.lat_lon = [row["LATITUDE"],row["LONGITUDE"]]
            output.append(out)
        
        self.stations = output

    def load_PWQMN(self):
        """
        Load the entire PWQMN database to the river class
        """
        import sqlite3
        conn = sqlite3.connect('PWQMN_SQLite2.db')
        sql = 'SELECT * FROM PWQMN_database'
        df = pd.read_sql(sql, conn,parse_dates='DATE')
        df['STN_DT'] = df['STATION'] + df['DATE'].dt.strftime(' %Y-%m-%d')
        self.database = df
    ## END OF CLASS river



















    