#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 16 20:24:26 2021

@author: sandy
"""

import numpy as np
import nmrglue as ng
import nmrglue.fileio.pipe as pipe
import nmrglue.process.pipe_proc as p
from urllib.request import urlopen as uReq
from bs4 import BeautifulSoup as soup
import requests
from selenium import webdriver
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.chrome.options import Options
from lxml import html
from selenium.webdriver.support.ui import Select
import pandas as pd
from sklearn import metrics
import re



class bond:
    def __init__(self,bondID,connection,bondOrder=1):
        self.bondID=bondID
        self.connection=connection
        self.bondOrder=bondOrder
        
    def __repr__(self):
        return "bondID:%s connection:%s bondOrder:%s"%(self.bondID, self.connection, self.bondOrder)
    

class atom:
    def __init__(self, atomID, connections, neighboors, atomType="C"):
        self.atomID=atomID
        self.atomType=atomType
        self.connections=connections
        self.neighboors=neighboors
        
    def __repr__(self):
        return "id:%s connection:%s neighboor:%s atomType:%s"%(self.atomID, self.connections, self.neighboors, self.atomType)
            
    def getTotalOrder(self):
        bondCounts=0
        for c in self.connections:
            bondCounts=bondCounts+int(c.bondOrder)
        return bondCounts
    
def calcH(structure,atom):
    if atom.atomType == 'C':
        Hcount=4-atom.getTotalOrder()
        #for n in atom.neighboors:
            #if getAtombyID(structure,n).atomType =='H':
                #Hcount+=1
        return Hcount
    else:
        if 'H' in atom.atomType and len(atom.atomType)>2:
            return 1+int(atom.atomType[len(atom.atomType)-1])
        elif 'H' in atom.atomType:
            return 1
        elif atom.atomType =='N':
            Hcount=3-atom.getTotalOrder()
            for n in atom.neighboors:
                if getAtombyID(structure,n).atomType =='H':
                    Hcount+=1
            return Hcount
        elif atom.atomType =='O':
            Hcount=2-atom.getTotalOrder()
            for n in atom.neighboors:
                if getAtombyID(structure,n).atomType =='H':
                    Hcount+=1
            return Hcount
        elif atom.atomType =='S':
            Hcount=2-atom.getTotalOrder()
            for n in atom.neighboors:
                if getAtombyID(structure,n).atomType =='H':
                    Hcount+=1
            return Hcount
        else:
            return 0
        
def getAtombyID(structure,atomID):
    for a in structure:
        if a.atomID==atomID:
            return a

def inputStructure():
    global strFile
    strFile=input("Enter the predict structure file path: ")

def inputSpectrum():
    print(40*"_")
    print(3*" "+"Select the value below to retrieve")
    print(40*"_")
    print('{:>23}'.format("H NMR[0]"))
    print('{:>23}'.format("C NMR[1]"))
    global idChoice
    idChoice = int(input("Enter a number of choice: "))
    if 0>idChoice>1:
        print(2*'\n'+38*'*')
        print("Incorrect Number Choice, Try Again.")
        print(38*'*'+2*'\n')
        inputSpectrum()
    global specFile
    specFile=input("Enter the product NMR spectrum file path: ")
    
def inputPeakInfo():
    print(40*"_")
    print(3*" "+"Select the value below to retrieve")
    print(40*"_")
    print('{:>23}'.format("H NMR[0]"))
    print('{:>23}'.format("C NMR[1]"))
    global idChoice
    idChoice = int(input("Enter a number of choice: "))
    if 0>idChoice>1:
        print(2*'\n'+38*'*')
        print("Incorrect Number Choice, Try Again.")
        print(38*'*'+2*'\n')
        inputSpectrum()
    print("Enter the peaks location of your spectrum in [] and sepreate them with ,")
    peakInfo=input("If the intergration appears to be 2, Please enter the peak two times.")
    return peakInfo
    
def strInfo(strFile):
    info=[]
    foundBound=False
    bondCount=-1
    bondInfo=[]
    
    foundAtom=False
    atomCount=-1
    atomInfo=[]

    
    with open(strFile,'r') as f:
        info=f.readlines()
        for i in info:
            if "/n><b"in i:
                foundBound=True
                bondCount=bondCount+1
                bondInfo.append([])
                foundAtom=True
            if foundBound:
                if "id=" in i:
                    bondInfo[bondCount].append(i.split('"')[1])
                elif "B=" in i:
                    bondInfo[bondCount].append(i.split('"')[1])
                elif "E=" in i:
                    bondInfo[bondCount].append(i.split('"')[1])
                elif "Order=" in i:
                    bondInfo[bondCount].append(i.split('"')[1])
            if foundAtom:
                if "id=" in i:
                    atomInfo[atomCount].append(i.split('"')[1])
                if "</s>" in i:
                    str1=i.split(">")
                    for a in str1:
                        if "/s" in a:
                            atomInfo[atomCount].append(a.split('<')[0])
            if "<b"in i:
                foundBound=True
                bondCount=bondCount+1
                bondInfo.append([])
                foundAtom=False
            elif "<n" in i:
                foundAtom=True
                atomCount=atomCount+1
                atomInfo.append([])
                foundBound=False
            elif "/>" in i:
                foundBound=False
                foundAtom=False
         
    bonds=[]
    for b in bondInfo:
        if len(b)==3:
            bonds.append(bond(b[0],[b[1],b[2]]))
        if len(b)==4:
            bonds.append(bond(b[0],[b[1],b[2]],b[3]))
            
    neighboorList, connectionList = bondsAnalysis(bonds)
    atomList=[]
    for i in range(len(atomInfo)):
        if len(atomInfo[i])==1:
            atomList.append(atom(atomInfo[i][0],connectionList[atomInfo[i][0]],
                                 neighboorList[atomInfo[i][0]]))
        elif len(atomInfo[i])==2:
            atomList.append(atom(atomInfo[i][0],connectionList[atomInfo[i][0]],
                                 neighboorList[atomInfo[i][0]],atomInfo[i][1]))
    #print(atomList)
    #print(makeMatrix(atomList))
    return makeMatrix(atomList)

    
def bondsAnalysis(bondList):
    neighboorList={}
    connectionList={}
    for b in bondList:
        for i in range(0,2):
            if neighboorList.get(b.connection[i],-1)==-1:
                neighboorList[b.connection[i]]=[b.connection[i-1]]
                connectionList[b.connection[i]]=[b]
            else:
                neighboorList[b.connection[i]].append(b.connection[i-1])
                connectionList[b.connection[i]].append(b)

    return neighboorList, connectionList

def makeMatrix (atomList):
    AdjM = np.zeros((len(atomList),len(atomList)))
    featureM = np.zeros((len(atomList),8))
    for a in range(len(atomList)):
        for i in range(len(atomList[a].neighboors)):
            for n in atomList:
                if not atomList[a] == n:
                    if n.atomID==atomList[a].neighboors[i]:
                        AdjM[a][atomList.index(n)]=atomList[a].connections[i].bondOrder
                        #print("AAA")
        
        if atomList[a].atomType == 'C':
            featureM[a][0]=1
            #print("BBB")
        elif 'N' in atomList[a].atomType:
            featureM[a][1]=1
        elif 'O' in atomList[a].atomType:
            featureM[a][2]=1
        elif 'S' in atomList[a].atomType:
            featureM[a][3]=1
        
        if calcH(atomList,atomList[a])==0:
            featureM[a][4]=1
        elif calcH(atomList,atomList[a])==1:
            featureM[a][5]=1
        elif calcH(atomList,atomList[a])==2:
            featureM[a][6]=1
            #print("CCC")
        elif calcH(atomList,atomList[a])==3:
            featureM[a][7]=1
    #return np.concatenate((np.array(AdjM,dtype="object"),np.array(featureM,dtype="object")),axis=1)  
    return np.concatenate((np.array(AdjM,dtype="object"), np.array(featureM,dtype="object")),axis=None)

def getDatabase(expression):
    option=webdriver.ChromeOptions()
    option.add_argument('headless')
    driver=webdriver.Chrome(options=option,executable_path='/Users/sandy/Downloads/chromedriver')
    driver.get('https://nmrshiftdb.nmr.uni-koeln.de/portal/media-type/html/user/anon/page/default.psml/js_pane/P-Search;jsessionid=8F4D9E5B1A5D93FFEEB7A0E50AA8CDFB')
    searchbox=driver.find_element_by_xpath('/html/body/table[2]/tbody/tr/td/table[3]/tbody/tr/td/table/tbody/tr/td/table/tbody/tr/td[1]/table/tbody/tr[2]/td/div/table/tbody/tr[2]/td[2]/form/input[1]')
    searchbox.send_keys(expression)
    searchButton = driver.find_element_by_xpath('/html/body/table[2]/tbody/tr/td/table[3]/tbody/tr/td/table/tbody/tr/td/table/tbody/tr/td[1]/table/tbody/tr[2]/td/div/table/tbody/tr[2]/td[2]/form/input[2]')   
    searchButton.click()
    my_url=driver.current_url
    #driver.close()
    
    uClient=uReq(my_url)
    page_html=uClient.read()
    uClient.close()
    page_soup=soup(page_html,"html.parser")
    pages=page_soup.findAll("td",{"class":"ContentStyleClass"})
    restultNum=int(pages[0].get_text().split("\n")[18].split(" ")[1])
    if restultNum % 10 == 0:
        pagenumber=restultNum//10
    else:
        pagenumber=restultNum//10 + 1
    resultPage=[my_url]
    if restultNum>10:
        i=0
        while i < pagenumber-1:
            a=14+i
            resultPage.append('https://nmrshiftdb.nmr.uni-koeln.de/'+page_soup.findAll("a")[a]["href"])
            i+=1
    databaseCInput=[]
    databaseHInput=[]
    databaseCOutput=[]
    CSpecNum=0
    databaseHOutput=[]
    HSpecNum=0
    for j in resultPage:
        #print(j)
        uClient=uReq(j)
        page_html=uClient.read()
        uClient.close()
        page_soup=soup(page_html,"html.parser")
        tables=page_soup.findAll("table")
        links=tables[16].findAll("a")
        for k in range(1, len(links)-1):
            driver.get('https://nmrshiftdb.nmr.uni-koeln.de/'+links[k]["href"])
            try:
                NMRtype=driver.find_element_by_css_selector('body > table:nth-child(2) > tbody > tr > td > table:nth-child(3) > tbody > tr > td > table > tbody > tr > td > table > tbody > tr > td:nth-child(2) > table > tbody > tr:nth-child(2) > td > div > table > tbody > tr:nth-child(2) > td.ContentStyleClass > table:nth-child(7) > tbody > tr:nth-child(1) > td:nth-child(1)').text
                try:
                    myclick=driver.find_element_by_css_selector('body > table:nth-child(2) > tbody > tr > td > table:nth-child(3) > tbody > tr > td > table > tbody > tr > td > table > tbody > tr > td:nth-child(2) > table > tbody > tr:nth-child(2) > td > div > table > tbody > tr:nth-child(2) > td.ContentStyleClass > table:nth-child(5) > tbody > tr > td > table:nth-child(1) > tbody > tr > td:nth-child(5) > a')
                    myclick.click()
                    mySelect=Select(driver.find_element_by_name("format"))
                    my = mySelect.select_by_visible_text("cml code")
                    myclick=driver.find_element_by_css_selector('body > table:nth-child(2) > tbody > tr > td > table:nth-child(3) > tbody > tr > td > table > tbody > tr > td > table > tbody > tr > td:nth-child(2) > table > tbody > tr:nth-child(2) > td > div > table > tbody > tr:nth-child(2) > td.ContentStyleClass > form:nth-child(9) > input[type=submit]')
                except:
                    myclick=driver.find_element_by_css_selector('body > table:nth-child(2) > tbody > tr > td > table:nth-child(3) > tbody > tr > td > table > tbody > tr > td > table > tbody > tr > td:nth-child(2) > table > tbody > tr:nth-child(2) > td > div > table > tbody > tr:nth-child(2) > td.ContentStyleClass > table:nth-child(7) > tbody > tr > td > table:nth-child(1) > tbody > tr > td:nth-child(5) > a')
                    myclick.click()
                    mySelect=Select(driver.find_element_by_name("format"))
                    my = mySelect.select_by_visible_text("cml code")
                    myclick=driver.find_element_by_css_selector('body > table:nth-child(2) > tbody > tr > td > table:nth-child(3) > tbody > tr > td > table > tbody > tr > td > table > tbody > tr > td:nth-child(2) > table > tbody > tr:nth-child(2) > td > div > table > tbody > tr:nth-child(2) > td.ContentStyleClass > form:nth-child(11) > input[type=submit]')
                myclick.click()
            except:
                try:
                    NMRtype=driver.find_element_by_css_selector('body > table:nth-child(2) > tbody > tr > td > table:nth-child(3) > tbody > tr > td > table > tbody > tr > td > table > tbody > tr > td:nth-child(2) > table > tbody > tr:nth-child(2) > td > div > table > tbody > tr:nth-child(2) > td.ContentStyleClass > table:nth-child(9) > tbody > tr:nth-child(1) > td:nth-child(1)').text + driver.find_element_by_css_selector('body > table:nth-child(2) > tbody > tr > td > table:nth-child(3) > tbody > tr > td > table > tbody > tr > td > table > tbody > tr > td:nth-child(2) > table > tbody > tr:nth-child(2) > td > div > table > tbody > tr:nth-child(2) > td.ContentStyleClass > table:nth-child(14) > tbody > tr:nth-child(1) > td:nth-child(1)').text
                    myclick=driver.find_element_by_css_selector('body > table:nth-child(2) > tbody > tr > td > table:nth-child(3) > tbody > tr > td > table > tbody > tr > td > table > tbody > tr > td:nth-child(2) > table > tbody > tr:nth-child(2) > td > div > table > tbody > tr:nth-child(2) > td.ContentStyleClass > table:nth-child(7) > tbody > tr > td > table:nth-child(1) > tbody > tr > td:nth-child(5) > a')
                    myclick.click()
                    mySelect=Select(driver.find_element_by_name("format"))
                    my= mySelect.select_by_visible_text("cml code")
                    myclick=driver.find_element_by_css_selector('body > table:nth-child(2) > tbody > tr > td > table:nth-child(3) > tbody > tr > td > table > tbody > tr > td > table > tbody > tr > td:nth-child(2) > table > tbody > tr:nth-child(2) > td > div > table > tbody > tr:nth-child(2) > td.ContentStyleClass > form:nth-child(11) > input[type=submit]')
                    myclick.click()
                except:
                    NMRtype=driver.find_element_by_css_selector('body > table:nth-child(2) > tbody > tr > td > table:nth-child(3) > tbody > tr > td > table > tbody > tr > td > table > tbody > tr > td:nth-child(2) > table > tbody > tr:nth-child(2) > td > div > table > tbody > tr:nth-child(2) > td.ContentStyleClass > table:nth-child(10) > tbody > tr:nth-child(1) > td:nth-child(1)').text
                    myclick=driver.find_element_by_css_selector('body > table:nth-child(2) > tbody > tr > td > table:nth-child(3) > tbody > tr > td > table > tbody > tr > td > table > tbody > tr > td:nth-child(2) > table > tbody > tr:nth-child(2) > td > div > table > tbody > tr:nth-child(2) > td.ContentStyleClass > table:nth-child(8) > tbody > tr > td > table:nth-child(1) > tbody > tr > td:nth-child(5) > a')
                    myclick.click()
                    mySelect=Select(driver.find_element_by_name("format"))
                    my = mySelect.select_by_visible_text("cml code")
                    myclick=driver.find_element_by_css_selector('body > table:nth-child(2) > tbody > tr > td > table:nth-child(3) > tbody > tr > td > table > tbody > tr > td > table > tbody > tr > td:nth-child(2) > table > tbody > tr:nth-child(2) > td > div > table > tbody > tr:nth-child(2) > td.ContentStyleClass > form:nth-child(12) > input[type=submit]')
                    myclick.click()
            driver.switch_to_window(driver.window_handles[-1])
            #print(NMRtype)
            if '1H' in NMRtype and '13C' in NMRtype:
                stucture,inputMatrix=getStructure(driver.current_url)
                typesList=list(NMRtype)
                if typesList.index('C')>typesList.index('H'):
                    Horder,Corder=1,2
                else:
                    Horder,Corder=2,1
                usable,NMRdata= getHNMRData(Horder,stucture,'https://nmrshiftdb.nmr.uni-koeln.de/'+links[k]["href"])
                if usable:
                    databaseHInput.append(inputMatrix)
                    databaseHOutput.append(NMRdata)
                    HSpecNum+=1
                usable,NMRdata= getCNMRData(Corder,stucture,'https://nmrshiftdb.nmr.uni-koeln.de/'+links[k]["href"])
                if usable:
                    databaseCInput.append(inputMatrix)
                    databaseCOutput.append(NMRdata)
                    CSpecNum+=1
            elif '1H' in NMRtype:
                stucture,inputMatrix=getStructure(driver.current_url)
                usable,NMRdata= getHNMRData(1,stucture,'https://nmrshiftdb.nmr.uni-koeln.de/'+links[k]["href"])
                if usable:
                    databaseHInput.append(inputMatrix)
                    databaseHOutput.append(NMRdata)
                    HSpecNum+=1
            elif '13C' in NMRtype:
                stucture,inputMatrix=getStructure(driver.current_url)
                usable,NMRdata= getCNMRData(1,stucture,'https://nmrshiftdb.nmr.uni-koeln.de/'+links[k]["href"])
                if usable:
                    #print(k)
                    databaseCInput.append(inputMatrix)
                    databaseCOutput.append(NMRdata)
                    #print(databaseCInput[CSpecNum].shape)
                    CSpecNum+=1
            #driver.close()
    #print(np.array(databaseHInput))
    #print(np.array(databaseHOutput))
    #print(np.reshape(np.array(databaseCInput,dtype="object"),[COutputNum*databaseCInput[0].shape[0],databaseCInput[0].shape[1]]))
    #print(np.array(databaseCOutput))
    #return np.array(databaseInput), np.array(databaseOutput)
    databaseCInput=np.array(databaseCInput,dtype="object")
    #print("nparray",type(databaseCInput))
    databaseHInput=np.array(databaseHInput,dtype="object")
    databaseCOutput=np.array(databaseCOutput,dtype="object")
    databaseHOutput=np.array(databaseHOutput,dtype="object")
    #print('There are ',CSpecNum,' C Spectrums')
    #print(databaseCInput)
    return [databaseHInput,databaseHOutput], [databaseCInput,databaseCOutput]
    #getNMRData(my_url)

def getNumC(atomList):
    numC=0
    for a in atomList:
        if a.atomType == 'C':
            numC+=1
    return numC

def getCNMRData(tableOrder,structure,url):
    uClient=uReq(url)
    page_html=uClient.read()
    uClient.close()
    page_soup=soup(page_html,"html.parser")
    if tableOrder ==1:
        table=page_soup.find("table",{"onmouseout":"select(-1,document.Spectrum1,document.JcpViewer1)"})
    else:
        table=page_soup.find("table",{"onmouseout":"select(-1,document.Spectrum2,document.JcpViewer2)"})
    trcontainer=table.findAll("tr")
    NMRshift=[]
    for tr in trcontainer:
        if(len(tr.findAll("td"))!= 0):
            NMRshift.append(float(tr.findAll("td")[2].get_text()))
    usable=True
    if len(NMRshift) < getNumC(structure):
        usable=False
    return usable, np.array(NMRshift,dtype="object")

def getHNMRData(tableOrder,structure,url):
    uClient=uReq(url)
    page_html=uClient.read()
    uClient.close()
    page_soup=soup(page_html,"html.parser")
    if tableOrder ==1:
        table=page_soup.find("table",{"onmouseout":"select(-1,document.Spectrum1,document.JcpViewer1)"})
    else:
        table=page_soup.find("table",{"onmouseout":"select(-1,document.Spectrum2,document.JcpViewer2)"})
    #print(table)
    trcontainer=table.findAll("tr")
    NMRshift=[]
    usable=True
    for a in structure:
        found=False
        for tr in trcontainer:
            if len(tr.findAll("td"))!=0:
                #print(tr.findAll("td")[0].get_text().replace(" ","").split("-")[0])
                #print(a.atomID)
                if int(tr.findAll("td")[0].get_text().replace(" ","").split("-")[0]) == int(a.atomID):
                    found=True
                    #print(a.atomID)
                    #print("Enter1")
                    if tr.findAll("td")[2].get_text()!= '':
                        data=tr.findAll("td")[2].get_text().split(" (")[0]
                        if ", " in data:
                            data=data.split(", ")
                            for d in data:
                                NMRshift.append(float(d))
                        else:
                            NMRshift.extend([float(data),0.0])
                    else:
                        if calcH(structure,a) ==0:
                            NMRshift.append(0.0)
                        else:
                            #print('Error 1')
                            usable=False
        if not found:
            if a.atomType !='H':
                if calcH(structure,a) ==0:
                    #print("Enter2")
                    NMRshift.append(0.0)
                    NMRshift.append(0.0)
                else:
                    #print('Error 2')
                    usable=False
    return usable, np.array(NMRshift,dtype="object")

def getStructure(url):
    r = requests.get(url, stream=True)
    bondList=[]
    atomInfo=[]
    for line in r.iter_lines():
        temp=str(line)
        if "bond " in temp:
            connection=[]
            for atomID in temp.split('"')[3].split(" "):
                connection.append("".join(filter(str.isdigit,atomID)))
            if temp.split('"')[5]=="S":
                bondList.append(bond(temp.split('"')[1],connection, 1))
            elif temp.split('"')[5]=="D":
                bondList.append(bond(temp.split('"')[1], connection, 2))
            elif temp.split('"')[5]=="T":
                bondList.append(bond(temp.split('"')[1], connection, 3))
        elif "atom " in temp:
            atomInfo.append(["".join(filter(str.isdigit,temp.split('"')[1])),temp.split('"')[3]])
    neighboorList, connectionList = bondsAnalysis(bondList)
    atomList=[]
    for i in range(len(atomInfo)):
        if(atomInfo[i][1]!="H"):
            atomList.append(atom(atomInfo[i][0],connectionList[atomInfo[i][0]],neighboorList[atomInfo[i][0]],atomInfo[i][1]))
    #return atomList
    return atomList, makeMatrix(atomList)

def readNMR():
    dic,data = ng.pipe.read(specFile)
    print (data[0])
    
class DeepNeuralNetwork:
    def __init__(self, layerIDToNumUnits, learning_rate, num_epochs, print_cost):
        self.layerIDToNumUnits = layerIDToNumUnits 
        self.L = len(self.layerIDToNumUnits)-1   
        self.nL = self.layerIDToNumUnits[self.L] 
        self.learning_rate = learning_rate
        self.num_epochs = num_epochs
        self.print_cost = print_cost
        self.costs = []
        
    def __initParams(self):
        self.params = {}
        rng = np.random.default_rng(42)
        for i in range(1, self.L+1):
            multiplier = np.sqrt(2/self.layerIDToNumUnits[i-1])
            
            self.params["W" + str(i)] = rng.standard_normal((self.layerIDToNumUnits[i], self.layerIDToNumUnits[i-1])) * multiplier # Note the required shapes!
            self.params["b" + str(i)] = np.zeros((self.layerIDToNumUnits[i], 1))
            
    def __updateParams(self):
        for layerID in range(1, self.L+1):
            sL = str(layerID)
            self.params["W"+sL] = self.params["W"+sL] - self.learning_rate * self.params["dW"+sL]
            self.params["b"+sL] = self.params["b"+sL] - self.learning_rate * self.params["db"+sL]
            
    def reportParams(self, maxDim=5):
        print("##### reportParams #####")
        for k, v in self.params.items():
            print("----------")
            print(k, v.shape)
            
            rBound = min(maxDim, v.shape[0])
            for r in range(rBound):
                
                cBound = min(maxDim, v.shape[1])
                for c in range(cBound):
                    print("{:12.8f}".format(v[r,c]), "  ", end="")
                    
                if (cBound < v.shape[1]):
                    print("...")
                else:
                    print()
                    
            if (rBound < v.shape[0]):
                print("...\n")
            else:
                print()
        print("##### end of reportParams #####")
        
    def __assertShapes(self):
        for layerID in range(1, self.L+1):
            sL = str(layerID)
            self.__assertOneShape("W"+sL)
            self.__assertOneShape("b"+sL)
            
    def __assertOneShape(self, s):
        assert(self.params["d"+s].shape == self.params[s].shape)
        
    def __sigmoid(self, z):
        #print(np.exp(-z))
        epsilon = 1.0e-12
        result = 1/(1+np.exp(-z) + epsilon) + epsilon
        #print(type(z))
        return result
    
    def __sigmoidPrime(self, Z):
        #print("this is Z",Z)
        return self.__sigmoid(Z)*(1-self.__sigmoid(Z))
    
    def __relu(self, Z):
        RL = []
        for i in Z:
            rl = []
            for z in i:
                rl.append(max(0,z))
            RL.append(rl)
        return np.array(RL)
    
    def __reluPrime(self, Z):
        RL = []
        for i in Z:
            rl = []
            for z in i:
                if (z<0):
                    dg=0.0
                else:
                    dg=1.0
                rl.append(dg)
            RL.append(rl)
        return np.array(RL)
    
    def __forwardOneLayer(self, layerID, activationFunc):
        self.params["Z"+str(layerID)] = np.dot(self.params["W"+str(layerID)],self.params["A"+str(layerID-1)])+self.params["b"+str(layerID)]
        self.params["A"+str(layerID)] = activationFunc(self.params["Z"+str(layerID)])
        
    def forward(self, X):
        self.params["A"+str(0)]=X
        for i in range(1,self.L+1):
            if(i<self.L):
                self.__forwardOneLayer(i,self.__relu)
            else:
                self.__forwardOneLayer(i,self.__sigmoid)
        return self.params["A"+str(self.L)]
    
    def __backwardOneLayer(self, layerID, activationFuncPrime, m):
        self.params["dZ"+str(layerID)] = self.params["dA"+str(layerID)]*activationFuncPrime(self.params["Z"+str(layerID)])
        self.params["dW"+str(layerID)] = np.dot(self.params["dZ"+str(layerID)],self.params["A"+str(layerID-1)].transpose())/m
        self.params["db"+str(layerID)] = np.sum(self.params["dZ"+str(layerID)],axis=1,keepdims=True)/m
        self.params["dA"+str(layerID-1)] = np.dot(self.params["W"+str(layerID)].transpose(),self.params["dZ"+str(layerID)])
        
    def __backward(self, Y, m):
        self.params["dA"+str(self.L)] = -np.divide(Y,self.params["A"+str(self.L)])+np.divide(1-Y,1-self.params["A"+str(self.L)])
        for l in range(self.L, 0, -1):
            if(l==self.L):
                self.__backwardOneLayer(l,self.__sigmoidPrime,m)
            else:
                self.__backwardOneLayer(l, self.__reluPrime, m)
                
    def __trackCost(self, Y, networkOutputs, epochID):
        if epochID % 100 == 0: # only do this every 100th epoch
            cost = self.__computeCost(Y, networkOutputs)
            self.costs.append(cost) # record the cost
            if self.print_cost:
                print("Cost after epoch %i: %f" %(epochID, cost))
                
    def __computeCost(self, Y, networkOutputs):
        cost = -np.sum(Y*np.log(networkOutputs)+(1-Y)*np.log(1-networkOutputs))/networkOutputs.size
        return cost
    
    def fit(self, X, Y):
        self.__initParams()
        self.costs = []
        for i in range(0,self.num_epochs):
            #self.__initParams()
            self.forward(X)
            self.__trackCost(Y,self.params["A"+str(self.L)],i)
            self.__backward(Y,X.shape[1])
            self.__updateParams()
        #self.reportParams()
            #print(self.__computeCost(Y,self.params["A"+str(self.L)]))
            
    def predict(self, X):
        predictions = self.forward(X)
        #self.reportParams()
        return predictions
    
    def test():
        print("========================================")
        print("========================================")
        print("========================================")
        print("---------- test_fit ----------")
                                             #              layer l    layer 2    layer 3    layer 4=L
        alg = DeepNeuralNetwork([4, 1], # 4 inputs, 5 hidden1, 3 hidden2, 4 hidden3, 1 output
                               learning_rate=0.2, num_epochs=100, print_cost=True)
        
        X = np.array([[30.0,  0.1,  0.2],  # 3 examples, 4 input attributes each
                      [0.8,  40.0,  0.4], 
                      [0.9,  0.5,  50.0], 
                      [0.1,  0.7,  60.0]])
        print(X)
        Y = np.array([[0.8, 0.0, 0.08]])
        
        alg.fit(X, Y)
        #alg1=DeepNNClassifier([4, 5, 3, 4, 1], # 4 inputs, 5 hidden1, 3 hidden2, 4 hidden3, 1 output
                               #learning_rate=0.02, num_epochs=5, print_cost=True)
        #alg.__initParams()
        print(alg.predict(np.array([[0.2],[0.4],[.6],[.8]])))
        #alg.reportParams(5)
    
def compareResults(Y_predicted, Y_correct):
    return metrics.accuracy_score(np.squeeze(Y_predicted), np.squeeze(Y_correct))

def score(prediction, actual):
    if prediction.shape==actual.shape:
        return 1.0-np.mean(np.abs((np.sort(prediction) - np.sort(actual)) / np.sort(prediction)))
    elif prediction.shape>actual.shape:
        for i in range(0,prediction.shape[0]-actual.shape[0]):
            actual=np.append(actual,[0.0])
        return (1.0-np.mean(np.abs((np.sort(prediction) - np.sort(actual)) / np.sort(prediction))))*100
    else:
        print("Wrong Spectrum!")
        return 0
    
def main():
    inputstr=input("Enter the training formula: ")
    inputStructure()
    StructureInfo=[]
    for i in strInfo(strFile):
        StructureInfo.append(i+0.01)
    StructureInfo=np.array(StructureInfo)
    trainSet=getDatabase(inputstr)
    CNMR_exist=True
    HNMR_exist=True
    if len(trainSet[0][0])==0:
        HNMR_exist=False
    if len(trainSet[1][0])==0:
        CNMR_exist=False
    if HNMR_exist:
        HNMRX_train=[]
        for i in trainSet[0][0].transpose():
            HNMRX_trainColumn=[]
            for y in i:
                HNMRX_trainColumn.append(y+0.01)
            HNMRX_train.append(HNMRX_trainColumn)
        HNMRX_train=np.array(HNMRX_train).transpose()
        HNMRY_train=trainSet[0][1].transpose()
        HNMRNetwork = DeepNeuralNetwork([HNMRX_train.shape[0],HNMRY_train.shape[0]+1,HNMRY_train.shape[0]],learning_rate=0.2, num_epochs=4000, print_cost=True)
        HNMRNetwork.fit(HNMRX_train,HNMRY_train)
        HNMRPrediction=HNMRNetwork.predict(StructureInfo)
        print("The prediction of HNMR Spectrum is: ",HNMRPrediction) 
    else:
        print("Not enough HNMR data!")
    
    if CNMR_exist:
        CNMRX_train=[]
        for i in trainSet[1][0].transpose():
            CNMRX_trainColumn=[]
            for y in i:
                CNMRX_trainColumn.append(y+0.01)
            CNMRX_train.append(CNMRX_trainColumn)
        CNMRX_train=np.array(CNMRX_train)
        CNMRY_train=trainSet[1][1].transpose()
        CNMRNetwork = DeepNeuralNetwork([CNMRX_train.shape[0],CNMRY_train.shape[0]+1,CNMRY_train.shape[0]],learning_rate=0.2, num_epochs=4000, print_cost=True)
        CNMRNetwork.fit(CNMRX_train,CNMRY_train)
        CNMRPrediction=CNMRNetwork.predict(StructureInfo)
        print("The prediction of CNMR Spectrum is: ",CNMRPrediction) 
    else:
        print("Not enough CNMR data!")
    
    knownSpecInfo=np.array(inputPeakInfo())
    if idChoice ==0:
        print("The score of this spectrum is similar to this structure is: ",score(HNMRPrediction,knownSpecInfo))
    elif idChoice==1:
        print("The score of this spectrum is similar to this structure is: ",score(CNMRPrediction,knownSpecInfo))
    
def test():
    inputstr=input("Enter the training formula: ")
    trainSet=getDatabase(inputstr)
    CNMRX_train=[]
    print(trainSet[1][0].transpose().shape)
    for i in trainSet[1][0].transpose():
        CNMRX_trainColumn=[]
        for y in i:
            CNMRX_trainColumn.append(y+0.01)
        CNMRX_train.append(CNMRX_trainColumn)
    CNMRX_train=np.array(CNMRX_train)
    CNMRY_pretrain=trainSet[1][1].transpose()/1000
    #print(CNMRY_train)
    X = np.array([[0.7,  0.1,  0.2],  # 3 examples, 4 input attributes each
                      [0.8,  0.3,  0.4], 
                      [0.9,  0.5,  0.6], 
                      [0.1,  0.7,  0.8]])
    Y = np.array([[1.03, 1.0, 0.0],
                  [1.0,1.0,0.80],
                  [1.0,1.0,0.0],
                  [1.0,1.08,0.0],
                  [1.09,1.0,0.0],
                  [1.0,1.0,0.90]])
    
    Y1=np.array([[0.1285, 0.15259999999999999, 0.004],
                 [0.1285, 0.12490000000000001, 0.072],
                 [0.1285, 0.1343, 0.0648],
                 [0.1285, 0.1343, 0.0648],
                 [0.1285, 0.12490000000000001, 0.072],
                 [0.1285, 0.12340000000000001, 0.004]])
    
    Y2=[]
    for i in CNMRY_pretrain:
        CNMRY_trainColumn=[]
        for y in i:
            CNMRY_trainColumn.append(i)
#        Y2.append(CNMRY_trainColumn)
        Y2=np.array(Y2)
    
    #CNMRNetwork = DeepNeuralNetwork([CNMRX_train.shape[0],CNMRY_train.shape[0]+1,CNMRY_train.shape[0]],learning_rate=0.2, num_epochs=4000, print_cost=True)
    CNMRNetwork = DeepNeuralNetwork([CNMRX_train.shape[0],10,6],learning_rate=0.2, num_epochs=4000, print_cost=True)
    #CNMRNetwork.fit(CNMRX_train,CNMRY_train)
    CNMRNetwork.fit(CNMRX_train,Y1)
    #print(CNMRNetwork.predict(np.array([[.2],[.4],[.6],[.8]])))
    inputStructure()
    StructureInfo=[]
    for i in strInfo(strFile):
        StructureInfo.append(i+0.01)
    StructureInfo=np.array(StructureInfo).reshape(CNMRX_train.shape[0],1)
    CNMRPrediction=CNMRNetwork.predict(StructureInfo).transpose()*1000
    print("The prediction of CNMR Spectrum is: ",CNMRPrediction)
    #knownSpecInfo=np.array(inputPeakInfo())
    #print("The score of this spectrum is similar to this structure is: ",score(CNMRPrediction,knownSpecInfo))

def testScoring():
    print("========================================")
    print("========================================")
    print("========================================")
    print("---------- testScoring ----------")  
    predict=np.array([152.6,124.9,134.3,134.3,124.9,123.4])
    act1=np.array([152.6,124.9,134.3,134.3,124.9,123.4])
    act2=np.array([142.6,124.9,134.3,134.3,122.9,123.4])
    act3=np.array([152.6,124.9,134.3,134.3,124.9])
    act4=np.array([152.6,124.9,134.3,134.3,124.9,123.4,123.4])
    act5=np.array([175.78045141, 145.87497999, 170.73612238, 181.4917545,  145.72094628,179.24648466])
    
    print("This is the score of the same:",score(predict,act1))
    print("This is the score of one different:",score(predict,act2))
    print("This is the score of one less:",score(predict,act3))
    print("This is the score of one more:",score(predict,act4))
    print("This is a acture example: ",score(predict,act5))

    

