import copy
import json
import math
import re
from decimal import Decimal
import time
import matplotlib.pyplot as plt


def load_cocoa_json(path):
    cocoa_list = []
    with open(path, 'r') as cocoa_json_file:
        cocoa_json = json.load(cocoa_json_file)
        exposure_windows = cocoa_json['exposure_windows']
        for window in exposure_windows:
            window['DateMillisSinceEpoch'] = time.strftime('%Y-%m-%d %H:%M:%S',
                                                           time.localtime(window['DateMillisSinceEpoch'] / 1000))
            window = drop_log(window)
            cocoa_list.append(window)
    return cocoa_list


def drop_log(window: dict):
    del window['CalibrationConfidence']
    del window['Infectiousness']
    del window['ReportType']
    return window


def date_filter(cocoa_json, y, m, d):
    year = y
    month = m
    day = d
    date = f'{year}-{month}-{day} 09:00:00'
    new_cocoa_list = []
    for window in cocoa_json:
        if window['DateMillisSinceEpoch'] == date:
            new_cocoa_list.append(window)
    return new_cocoa_list


def dB2dis(cocoa_json):
    new_cocoa_json = copy.deepcopy(cocoa_json)
    for window in new_cocoa_json:
        scan_instances = window['ScanInstances']
        for instance in scan_instances:
            instance['MinAttenuationDb'] = f"{d2d(instance['MinAttenuationDb'])} cm"
            instance['TypicalAttenuationDb'] = f"{d2d(instance['TypicalAttenuationDb'])} cm"
    return new_cocoa_json


def d2d(db):
    try:
        db = int(db)
    except ValueError:
        db = float(db)
    distance = 0.0981 * math.pow(math.e, 0.1174 * db)
    return float(Decimal(distance).quantize(Decimal('0.1'), rounding='ROUND_HALF_UP'))


def summary_distance(cocoa_json):
    new_cocoa_json = copy.deepcopy(cocoa_json)
    for window in new_cocoa_json:
        total_time = 0
        total_db = 0
        min_db = 10000
        scan_instances = window['ScanInstances']
        for instance in scan_instances:
            total_time += instance['SecondsSinceLastScan']
            total_db += int(instance['SecondsSinceLastScan']) * int(instance['TypicalAttenuationDb'])
            if instance['MinAttenuationDb'] < min_db:
                min_db = instance['MinAttenuationDb']
        window['ScanInstances'] = {'MinAttenuationDistance': f'{d2d(min_db)} cm',
                                   'TotalTime': f'{total_time / 60} mins',
                                   'TypicalAttenuationDb': f'{d2d(total_db / total_time)} cm'}

    return new_cocoa_json


def draw_distance_log(cocoa_list):
    window_str = input('Paste a window line or type the No: \n').replace("'", '"')
    if window_str.isdigit():
        window_str = f'{cocoa_list[int(window_str)]}'.replace("'", '"')
        print(window_str)
    else:
        window_str = re.sub('No.\d+', '', window_str)
    window = json.loads(window_str)
    total_time = 0
    x_time = []
    y_min_distance = []
    y_avg_distance = []
    scan_instances = window['ScanInstances']
    for instance in scan_instances:
        total_time += instance['SecondsSinceLastScan'] / 60
        x_time.append(total_time)
        min_dis = 500
        avg_dis = 500
        if d2d(instance['MinAttenuationDb']) < min_dis:
            min_dis = d2d(instance['MinAttenuationDb'])
        if d2d(instance['TypicalAttenuationDb']) < avg_dis:
            avg_dis = d2d(instance['TypicalAttenuationDb'])
        y_min_distance.append(min_dis)
        y_avg_distance.append(avg_dis)
    plt.scatter(x_time, y_min_distance, c='r')
    plt.scatter(x_time, y_avg_distance, c='b')
    plt.plot(x_time, y_min_distance, 'r--', label='min distance')
    plt.plot(x_time, y_avg_distance, 'b--', label='avg distance')
    for x, y in zip(x_time, y_min_distance):
        plt.text(x, y, y)
    for x, y in zip(x_time, y_avg_distance):
        plt.text(x, y, y)
    plt.xlabel('time [min]')
    plt.ylabel('distence [cm]')
    plt.legend()
    plt.show()


if __name__ == '__main__':
    path = "X:\XXXXX\exposure_data.json"
    cocoa_list1 = load_cocoa_json(path)
    cocoa_list1 = date_filter(cocoa_list1, '2022', '07', '22')
    cocoa_list2 = dB2dis(cocoa_list1)
    cocoa_list3 = summary_distance(cocoa_list1)
    for i in range(len(cocoa_list1)):
        print(f'No.{i}', cocoa_list1[i])
        print(f'No.{i}', cocoa_list2[i])
        print(f'No.{i}', cocoa_list3[i])
        print('\n')
    print('\n\n')

    while True:
        draw_distance_log(cocoa_list1)
