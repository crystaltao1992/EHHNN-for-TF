import xlwt
import xlrd
import numpy

# filenames = ["311974","313740", "314559"]
# filenames = ["311974","313658", "314559"]

class DataSets(object):
    def std_data(self,data):
        for attr in data:
            max_flow = max(attr)
            min_flow = min(attr)
            gap = max_flow - min_flow
            # print(max_flow,min_flow,'max min')
            for i in range(len(attr)):
                attr[i] = (attr[i] - min_flow) / gap


    def savedata(self):
        outfilename = 'result/traffic_'+str(self.lag)+'_'+str(self.step)+'_'+str(self.station_num)+'.data'
        with open(outfilename, 'w') as f:
            for line in self.data:
                line = [str(each) for each in line]
                res = " ".join(line)
                f.write(res + '\n')


def get_col(station,attribute='Flow (Veh/5 Minutes)'):
    col_data=[]
    for i in range(-2,2):
        if i:
            filename= station+'_' +str(i)
        else:
            filename=station
        f_name = 'data//' + filename + '.xlsx'
        print(f_name)
        workbook = xlrd.open_workbook(f_name)
        booksheet = workbook.sheet_by_index(0)  # 用索引取第一个sheet
        row0 = booksheet.row_values(0)
        attribute_index = row0.index(attribute)
        col_data += booksheet.col_values(attribute_index)[1:]
    return col_data


def generate_traffic_data(lag=10, step=1, stations=[]):
    datasets = DataSets()
    datasets.step=step
    datasets.lag=lag
    datasets.station_num = len(stations)
    datasets.flows = []
    datasets.density=[]
    datasets.speed = []
    for each in stations:
        datasets.flows.append(get_col(each))
        datasets.density.append(get_col(each,'Occupancy (%)'))
        datasets.speed.append(get_col(each, 'Speed (mph)'))

    datasets.std_data(datasets.flows)
    datasets.std_data(datasets.density)
    datasets.std_data(datasets.speed)
    datasets.len = len(datasets.flows[0])
    datasets.data = []
    for i in range(datasets.len - lag- step):
        data = []
        station_num = len(datasets.flows)
        for station_id in range(station_num):
            data += datasets.speed[station_id][i:i + lag ]
            data += datasets.density[station_id][i:i + lag]
            data += datasets.flows[station_id][i :i + lag]
            if station_id==station_num-1:
                data += datasets.flows[station_id][i+lag +step-1:i+lag+step]
        datasets.data.append(data)
    return datasets


if __name__ == '__main__':
    lag, step = 10, 1
    stations = ['311974', '318282', '314270', '319631', '313852', '313658', '314559']
    datasets = generate_traffic_data(lag=lag,step=step,stations=stations)
    datasets.savedata()
