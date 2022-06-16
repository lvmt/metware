#!/usr/bin/python3
import re, sys
import lxml
import requests
from bs4 import BeautifulSoup
pyfile, outfile = sys.argv

class kegg(object):
    def __init__(self):
        self.url = 'https://www.genome.jp/kegg/pathway.html'
        self.headers = {'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/87.0.4280.67 Safari/537.36 Edg/87.0.664.55'}

    def get_html(self, url):
        response = requests.get(url, headers=self.headers)
        if response.status_code == 200:
            return response
        else:
            return None

    def keggPathway(self, url):
        if self.get_html(url):
            html = self.get_html(url).text
            soup = BeautifulSoup(html, "lxml")
            Kmain = soup.find(class_="main")
            # print(Kmain)
            fp = []
            sp = []
            for i in Kmain:
                if re.findall(r"<h4 id=.*>(.*?)</h4>", str(i)):
                    fp.append(re.findall(r"<h4 id=.*>(.*?)</h4>", str(i)))
                if re.findall(r"<b.*>(.*?)</b>", str(i)):
                    sp.append(re.findall(r"<b.*>(.*?)</b>", str(i)))
            del sp[0]
            # print(fp,sp)
            kt = Kmain.find_all(class_="list")
            print(len(fp), len(sp), len(kt))
            tpw = []
            mapc = 0
            for i in fp:
                # print(i[0])
                for j, k in zip(sp, kt):
                    # print(j[0])
                    if i[0][0] == j[0][0]:
                        fs = re.sub(r'^[^A-Z]*', '', i[0]) + "\t" + re.sub(
                            r'^[^A-Z]*', '', j[0])
                        ka = k.find_all('a')

                        for l in ka:
                            kau = l.get("href")
                            kat = l.get_text()
                            mapN = re.findall(r'(\w\w\w=?\d{5})', str(kau))[0]
                            mapN = re.sub(r'=', '', mapN)
                            if "map" in mapN:
                                mapc += 1
                            tpw.append(fs + "\t" + kat + "\t" + mapN)
            print("map has :", mapc)
            return tpw

    def main(self):
        tpw = self.keggPathway(self.url)
        with open(outfile, "w") as f:
            f.write("level1\tlevel2\tlevel3\tmap\n")
            for i in tpw:
                f.write(i + "\n")

if __name__ == "__main__":
    kp = kegg()
    kp.main()
