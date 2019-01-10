import os
import re
import sys
import time
import json
import subprocess
import threading
from queue import Queue

import requests
import selenium
from lxml import etree
from requests.compat import urljoin
from selenium import webdriver
from selenium.webdriver.common.keys import Keys
from pyvirtualdisplay import Display




class UrlPool:
    """获取给定项目号的所有item，并保存到队列"""
    def __init__(self, queue, accession_num, startURL):
        self.accession_num = accession_num
        self.startURL = startURL
        self.queue = queue

    def parse_page(self, content):
        html = etree.HTML(content)
        node_list = html.xpath("//div[@class='rslt']/p[@class='title']/a/@href")
        if node_list:
            for node in node_list:
                self.queue.put(urljoin(self.startURL, node))

    def browse_page(self):
        # 启动虚拟桌面
        display = Display(visible=0, size=[800, 600])
        display.start()
        driver = webdriver.Chrome()
        driver.get(self.startURL)
        elem = driver.find_element_by_xpath('//input[@title="Search SRA"]')
        elem.clear()
        try:
            print(f'正在搜索{self.accession_num}')
            elem.send_keys(self.accession_num)
            elem.send_keys(Keys.RETURN)
        except Exception:
            raise('SRA检索失败，请检查项目号后再次运行。')

        has_next = True

        while has_next:
            try:
                # 将页面源代码提供给parse_page进行解析
                page_content = driver.page_source
                self.parse_page(page_content)
                # 点击<下一页>按钮进行翻页
                options = driver.find_elements_by_xpath("//a[@title='Next page of results']")
                options[0].click()
            except IndexError:
                has_next = False
        driver.close()
        display.stop()


class ThreadCrawl(threading.Thread):
    def __init__(self, thread_name, data_queue, output_dir):
        super(ThreadCrawl, self).__init__()
        self.thread_name = thread_name
        self.data_queue = data_queue
        self.output_dir = output_dir

    def run(self):
        print(self.thread_name, '正在作业...')
        while not PARSE_EXIT:
            try:
                url = self.data_queue.get(False)
                self.parse(url)
            except Exception:
                pass

    def parse(self, url):
        page = requests.get(url).content
        html = etree.HTML(page)
        srr = html.xpath("//td/a/text()")
        if srr:
            self.download(srr[0])

    def download(self, srr):
        print(f'{self.thread_name}正在下载{srr}')
        srrbase = re.search(r'SRR(\d{3})', srr).group()
        ascp_ = 'anonftp@ftp.ncbi.nlm.nih.gov:/sra/sra-instant/reads/ByRun/sra/SRR'
        try:
            ret = subprocess.run([
                    '/home/woods/.aspera/connect/bin/ascp', '-i',
                    '/home/woods/.aspera/connect/etc/asperaweb_id_dsa.openssh',
                    '-k','1','-T','-l200m',
                    os.path.join(ascp_, srrbase, srr, f'{srr}.sra'),
                    self.output_dir
                    ], stderr=subprocess.STDOUT)
            if ret.returcode != 0:
                print(f'Error downloading {srr},', ret)
        except OSError:
            pass

PARSE_EXIT = False

def main():
    accession_num = input('请输入项目号:')
    output_dir = input('请输入输出路径:')
    if not output_dir:
        output_dir = './'
    #output_dir = '/home/pub/HDisk/RawData/Microbiom/SRP059928'
    data_queue = Queue()
    url = 'https://www.ncbi.nlm.nih.gov/sra'
    print('默认输出路径为当前目录',output_dir)

    urlpool = UrlPool(data_queue, accession_num, url)
    urlpool.browse_page()

    crawlist = ['挖掘机1号', '挖掘机2号', '挖掘机3号', '挖掘机4号', '挖掘机5号']
    threadcrawl = []
    for crawl in crawlist:
        thread = ThreadCrawl(crawl, data_queue, output_dir)
        thread.start()
        threadcrawl.append(thread)

    while not data_queue.empty():
        pass

    global PARSE_EXIT
    PARSE_EXIT = True

    for thread in threadcrawl:
        thread.join()
    
    print("数据下载完毕.")


if __name__ == '__main__':
    main()

        

