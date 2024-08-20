def extract_data_from_line(line):
    # 将每一行按逗号分割
    parts = line.split(',')

    # 检查是否有足够的数据
    if len(parts) < 15:
        return None

    # 提取所需的数据
    week = parts[4]
    wn = str(int(parts[5]) // 1000)
    prn = str(int(parts[10]) - 160)
    iodssr = parts[11]
    iodp = parts[12]
    todbdt = parts[13]

    # 将第15个数据转换为二进制，并去除'*'后的部分，确保长度为原十六进制字符串长度的4倍
    hex_data = parts[14].split('*')[0]
    binary_data = bin(int(hex_data, 16))[2:].zfill(len(hex_data) * 4)

    # 切片提取不同卫星系统的数据
    bds = binary_data[:63]
    gps = binary_data[63:100]
    galileo = binary_data[100:137]
    glonass = binary_data[137:174]

    return f'week:{week} wn:{wn} prn:{prn} iodssr:{iodssr} iodp:{iodp} todbdt:{todbdt} BDS:{bds} GPS:{gps} Galileo:{galileo} GLONASS:{glonass}\n'

# 读取文件并逐行处理
with open('D:\\cssrlibccc\\src\\cssrlib\\msg1.txt', 'r') as file, open('D:\\cssrlibccc\\src\\cssrlib\\bdsmsg1.txt', 'w') as output_file:
    for line in file:
        data = extract_data_from_line(line)
        if data:
            output_file.write(data)



