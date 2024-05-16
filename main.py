# def data_reader(root):
#     f = open(root)
#     data = f.readlines()
#     f.close()
#     return data

import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import chirp, find_peaks, peak_widths
from scipy.ndimage import gaussian_filter1d
import os
import matplotlib
matplotlib.rc("font",family='Dengxian')

class HPLC_data:
    def __init__(self, root, name='default'):
        f = open(root)
        data = f.readlines()
        f.close()

        self.name = name

        self.time_table = []
        for i in range(2, len(data)):
            self.time_table.append(float(data[i].split('\t')[0]))

        self.lamda_table = []
        for i in range(1, len(data[0].replace(' ', '').replace('\n', '').split('\t'))):
            self.lamda_table.append(float(data[0].replace(' ', '').replace('\n', '').split('\t')[i]))

        self.abs_table = []
        for i in range(2, len(data)):
            self.abs_table.append(list(map(float, data[i].replace(' ', '').replace('\n', '').split('\t')[1:])))

        self.time_table = np.array(self.time_table)
        self.lamda_table = np.array(self.lamda_table)
        self.abs_table = np.array(self.abs_table).T  # 初始化数据为可读的。吸收列表的横坐标为波长，纵坐标为时间。[lamda,time]

    # def draw_lamda(self, lamda=245, method='linear'):
    #     if method == 'linear':
    #         lamda_No = 0
    #         for i in range(len(self.lamda_table)):
    #             if lamda < self.lamda_table[i]:
    #                 lamda_No = i - 1
    #                 break
    #         former_abs = self.abs_table[lamda_No]
    #         latter_abs = self.abs_table[lamda_No + 1]
    #         absorb_table = []
    #         for i in range(len(former_abs)):
    #             absorb = (former_abs[i] * (self.lamda_table[lamda_No + 1] - lamda) + latter_abs[i] * (
    #                     lamda - self.lamda_table[lamda_No])) / (
    #                              self.lamda_table[lamda_No + 1] - self.lamda_table[lamda_No])  # 插值确定吸收曲线
    #             absorb_table.append(absorb)
    #
    #         # absorb_table = gaussian_filter1d(absorb_table,sigma=0.2)  # 曲线的平滑
    #
    #         peaks, _ = find_peaks(absorb_table[:6000], prominence=0.05)  # 寻峰。prominence为突出程度，越高，则只有越突出的峰会被记录。
    #         peaks2, _ = find_peaks(absorb_table[6001:], prominence=0.2)  # 寻峰。prominence为突出程度，越高，则只有越突出的峰会被记录。
    #         peaks = peaks + peaks2
    #         print(peaks)
    #
    #         plt.plot(self.time_table, absorb_table)
    #
    #         plt.plot(self.time_table[peaks], [absorb_table[i] for i in peaks], 'x')
    #         for i in peaks:
    #             plt.text(self.time_table[i], absorb_table[i],
    #                      '({:.3f},{:.3f})'.format(self.time_table[i], absorb_table[i]))
    #
    #         plt.show()
    #
    # def draw_time(self, time, method='linear'):
    #     if method == 'linear':
    #         time_No = 0
    #         for i in range(len(self.time_table)):
    #             if time < self.time_table[i]:
    #                 time_No = i - 1
    #                 break
    #
    #         absorb_table = []
    #         for i in range(len(self.abs_table)):
    #             absorb = (self.abs_table[i][time_No] * (self.time_table[time_No + 1] - time) + self.abs_table[i][
    #                 time_No + 1] * (time - self.time_table[time_No])) / \
    #                      (self.time_table[time_No + 1] - self.time_table[time_No])  # 插值确定吸收曲线
    #             absorb_table.append(absorb)
    #
    #         print(absorb_table)
    #
    #         # absorb_table = gaussian_filter1d(absorb_table,sigma=0.2)  # 曲线的平滑
    #
    #         peaks, _ = find_peaks(absorb_table,prominence=0.05)  # 寻峰。prominence为突出程度，越高，则只有越突出的峰会被记录。
    #         print(peaks)
    #
    #         plt.plot(self.lamda_table, absorb_table)
    #
    #         plt.plot(self.lamda_table[peaks], [absorb_table[i] for i in peaks], 'x')
    #         for i in peaks:
    #             plt.text(self.lamda_table[i], absorb_table[i],
    #                      '({:.3f},{:.3f})'.format(self.lamda_table[i], absorb_table[i]))
    #
    #         plt.show()

    def auto_analy_lamba(self, lamda=245, method='linear'):
        fig = plt.figure()
        if method == 'linear':
            lamda_No = 0
            for i in range(len(self.lamda_table)):
                if lamda < self.lamda_table[i]:
                    lamda_No = i - 1
                    break
            former_abs = self.abs_table[lamda_No]
            latter_abs = self.abs_table[lamda_No + 1]
            absorb_table = []
            for i in range(len(former_abs)):
                absorb = (former_abs[i] * (self.lamda_table[lamda_No + 1] - lamda) + latter_abs[i] * (
                        lamda - self.lamda_table[lamda_No])) / (
                                 self.lamda_table[lamda_No + 1] - self.lamda_table[lamda_No])  # 插值确定吸收曲线
                absorb_table.append(absorb)

            # absorb_table = gaussian_filter1d(absorb_table,sigma=0.2)  # 曲线的平滑
            # peaks, _ = find_peaks(absorb_table, prominence=0.2)  # 寻峰。prominence为突出程度，越高，则只有越突出的峰会被记录。
            #
            peaks, _ = find_peaks(absorb_table[:4000], prominence=0.05)  # 寻峰。prominence为突出程度，越高，则只有越突出的峰会被记录。
            peaks2, _ = find_peaks(absorb_table[4001:], prominence=0.2)  # 寻峰。prominence为突出程度，越高，则只有越突出的峰会被记录。
            peaks = np.append(peaks, (peaks2 + 4001))
            print(peaks)

            plt.xlim(0, 25)
            plt.ylim(-0.05, 1.5)
            plt.plot(self.time_table, absorb_table)
            plt.plot(self.time_table[peaks], [absorb_table[i] for i in peaks], 'x', color='purple')
            plt.title('{}, {} = {}nm'.format(self.name, chr(955), lamda))
            plt.xlabel('时间 / 分钟')
            plt.ylabel('吸光度')

            left, bottom, width, height = 0.15, 0.75, 0.08, 0.08
            axs = []
            for i in peaks:
                if self.time_table[i] > 18:
                    xx = self.time_table[i] - 5.5
                else:
                    xx = self.time_table[i]
                if absorb_table[i] > 1.1:
                    yy = 1.1
                else:
                    yy = absorb_table[i]

                plt.text(xx, yy, '({:.3f},{:.3f})'.format(self.time_table[i], absorb_table[i]))
            for i in peaks:
                time = self.time_table[i]
                if method == 'linear':
                    time_No = 0
                    for i in range(len(self.time_table)):
                        if time < self.time_table[i]:
                            time_No = i - 1
                            break

                    absorb_table_time = []
                    for i in range(len(self.abs_table)):
                        absorb = (self.abs_table[i][time_No] * (self.time_table[time_No + 1] - time) +
                                  self.abs_table[i][
                                      time_No + 1] * (time - self.time_table[time_No])) / \
                                 (self.time_table[time_No + 1] - self.time_table[time_No])  # 插值确定吸收曲线
                        absorb_table_time.append(absorb)
                    # absorb_table = gaussian_filter1d(absorb_table,sigma=0.2)  # 曲线的平滑
                    peaks, _ = find_peaks(absorb_table_time, prominence=0.05,
                                          distance=15)  # 寻峰。prominence为突出程度，越高，则只有越突出的峰会被记录。

                    axs.append(fig.add_axes([left, bottom, width, height]))
                    left += 0.12
                    axs[-1].plot(self.lamda_table, absorb_table_time, 'orange')

                    axs[-1].plot(self.lamda_table[peaks], [absorb_table_time[i] for i in peaks], 'x', color='red')
                    axs[-1].set_title('{:.3f}'.format(time), y=-0.6)
                    for i in peaks:
                        axs[-1].text(self.lamda_table[i], absorb_table_time[i],
                                     '{:.3f}'.format(self.lamda_table[i]))
                    axs[-1].axis('off')

            plt.savefig('./image\\{}.jpg'.format(self.name), dpi=750)
            # plt.show()


#
name_list = os.listdir('data')
for name in name_list:
    HPLC = HPLC_data('data\\'+name,name=name.split('.')[0])
    HPLC.auto_analy_lamba(245)
# A5 = HPLC_data('data\\G-12.txt')
# A5.auto_analy_lamba()
