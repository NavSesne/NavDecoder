def extract_data_from_line(line):
    parts = line.split(',')
    if len(parts) < 56:  # 确保每行有足够的数据
        return None
    week=parts[4]
    wn = str(int(parts[5]) // 1000)
    prn = str(int(parts[10]) - 160)
    iodssr = parts[11]
    iodp = parts[12]
    todbdt = parts[13]
    result = f'week:{week} wn:{wn} prn:{prn} iodssr:{iodssr} iodp:{iodp} todbdt:{todbdt}'

    for i in range(14, len(parts) - 1, 7):
        satslot = parts[i]
        iodn = parts[i + 1]
        Rorb = float(parts[i + 2]) * 0.0016
        Aorb = float(parts[i + 3]) * 0.0064
        Corb = float(parts[i + 4]) * 0.0064
        iodcorr = parts[i + 5]
        URAI = parts[i + 6].split('*')[0]  # 去除可能存在的'*'符号
        result += f' satslot:{satslot} iodn:{iodn} Rorb:{Rorb:.4f} Aorb:{Aorb:.4f} Corb:{Corb:.4f} iodcorr:{iodcorr} URAI:{URAI}'

    return result + '\n'

with open('msg2.txt', 'r') as file, open('bdsmsg2.txt', 'w') as output_file:
    for line in file:
        data = extract_data_from_line(line)
        if data:
            output_file.write(data)




