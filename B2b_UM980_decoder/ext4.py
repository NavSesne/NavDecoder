def extract_data_from_line(line):
    parts = line.split(',')
    if len(parts) < 56:  # 确保每行有足够的数据
        return None
    week = parts[4]
    wn = str(int(parts[5]) // 1000)
    prn = str(int(parts[10]) - 160)
    iodssr = parts[11]
    iodp = parts[12]
    todbdt = parts[13]
    sub = parts[14]
    result = f'week:{week} wn:{wn} prn:{prn} sub:{sub} iodssr:{iodssr} iodp:{iodp} todbdt:{todbdt}'
    for i in range(18, len(parts) - 1, 2):
        iodcorr = parts[i]
        sc0 = float(parts[i + 1]) * 0.0016
        result += f' iodcorr:{iodcorr} sc0:{sc0}'

    return result + '\n'


with open('msg4.txt', 'r') as file, open('bdsmsg4.txt', 'w') as output_file:
    for line in file:
        data = extract_data_from_line(line)
        if data:
            output_file.write(data)
