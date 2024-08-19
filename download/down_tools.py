#! python3
"""
@Description 	:   PPP-RTK processing/plotting library for POS/TRO/IONO
@Author      	:   Lewen Zhao
@Version     	:   1.0
@Contact     	:   lwzhao@nuist.edu.cn
@Contact     	:   navsense_support@163.com
@OnlineService :	http://1.13.180.60:8800/login2
@Time        	:   2024/03/29 12:14:41
"""
from datetime import datetime
import os
from dateutil.relativedelta import relativedelta
currentpath=os.path.dirname(os.path.abspath(__file__))
root_path=os.path.dirname(currentpath)
bin_wget=os.path.join(root_path,"bin",'wget.exe')
bin_curl=os.path.join(root_path,"bin",'curl.exe')
bin_gzip=os.path.join(root_path,"bin",'gzip.exe')
bin_crx =os.path.join(root_path,"bin",'crx2rnx.exe')
bin_gfzrnx=os.path.join(root_path,"bin",'gfzrnx_x64.exe')

def ymd2doy(year, mon, day):
    dn = datetime(year, mon, day, 0, 0, 0)
    return int(dn.strftime("%j"))

def call_wget_(dir_dst, url):
    if not os.path.exists(dir_dst):
        os.makedirs(dir_dst)
    cmd = '%s  -P %s %s' % (bin_wget,dir_dst, url)
    print (cmd)
    os.system(cmd)

def call_curl_(dir_dst, url):
    if not os.path.exists(dir_dst):
        os.makedirs(dir_dst)
    os.chdir(dir_dst)
    cmd = '%s  -s -S -c .cookie -n -L -O %s' % (bin_curl,url)
    print ("curl cmd="+cmd)
    os.system(cmd)
def download_compress(rpath,local_dir,name,suffix):
    # 检查本地是否存在文件
    compressed_name = name + suffix
    local_file_path = os.path.join(local_dir, name)
    local_compressed_file_path = os.path.join(local_dir, compressed_name)
    
    if os.path.exists(local_file_path):
        print(f"File {name} already exists.")
        return local_file_path
    if os.path.exists(local_compressed_file_path):
        os.remove(local_compressed_file_path)
        print(f"Compressed file {compressed_name} found and deleted.")
    
    # call_wget_(dir_dst, obsfs)
    remote_path=f"{rpath}{compressed_name}"
    call_curl_(local_dir, remote_path)    
    if os.path.exists(local_compressed_file_path):
        cmd='%s -d %s' %(bin_gzip,local_compressed_file_path)
        os.system(cmd)
    return local_file_path