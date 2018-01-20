# -*- coding: utf-8 -*-
"""
Created on Mon Nov 13 19:02:29 2017

@author: Big Pigeon
Note: https://phylogeographer.com/scripts/anonymizer.py
"""

# -*- coding: utf-8 -*-
"""
Created on Tue Oct 17 21:59:34 2017

@author: Big Pigeon
"""

import pandas as pd
import xml.etree.ElementTree
import sys

codesFile = sys.argv[1]

hierarchy = {}

paternalAncestryXMLFile = sys.argv[2]

idRegionOverrides = []

def parseOverridesFile():
    with open(overridesFile) as f:
        content = f.read().splitlines()
        for i in range(len(content)):
            splits = content[i].split(",")
            idRegionOverrides.append({'id': splits[0], 'replacement': splits[1]})
    print (idRegionOverrides)
    
if len(sys.argv) == 4:
    overridesFile = sys.argv[3]
    parseOverridesFile()

def parseHierarchyNov():
    with open(codesFile) as f:
        content = f.read().splitlines()
        for i in range(len(content)):
            splits = content[i].split(",")
            parent = splits[0]
            for j in range(len(splits) - 1):
                hierarchy[splits[j + 1]]=parent
    print (hierarchy)
    
parseHierarchyNov()

def parseSubgroupText(text):
    splits = text.split("&gt;")    
    for index in range (1, len(splits) + 1):
        thisCode = splits[len(splits) - index]        
        thisCode = thisCode.replace(" Confirmed/Predicted", "")
        if thisCode in hierarchy:
            return thisCode        
        else:
            questionSplit = thisCode.split("?")
            if questionSplit[0] in hierarchy:
                if len(questionSplit) == 3:
                    print ("two question marks, ignoring both guesses for", text)
                    return thisCode
                else:
                    print ("taking guess", questionSplit[0], "for", text)
                    return questionSplit[0]
            else:
                negativeSplit = thisCode.split("(x")                
                if len(negativeSplit) > 1:                    
                    if negativeSplit[0] in hierarchy:
                        return negativeSplit[0]
                    else:
                        trimmed = negativeSplit[0][0:-1]
                        if trimmed in hierarchy:
                            return trimmed
    return None

def parseSubgroupNegatives(text):
    negativeSplit = text.split("(x")
    if len(negativeSplit) > 1:
        closeParen = negativeSplit[1].index(')')
        negatives = negativeSplit[1][0:closeParen].split(",")
        for negative in negatives:
            if negative not in hierarchy:
                print("negative", negative, "not recognized")
        return negatives
    return []
    
parseReport = {}

    
def getValidKits():
        
    unknownLocations = set([])

    root = xml.etree.ElementTree.parse(paternalAncestryXMLFile).getroot()
    kitRows = root[1][0]
    validKits = []
                       
    countryCenters = {"Poland": [52.069, 19.480],
                      "Germany": [51.159, 10.448],
                      "Spain": [40.432, -3.703],
                      "England": [52.861, -1.223],
                      "Scotland": [56.034, -3.897],
                      "United Kingdom": [52.861, -1.223],
                      "Russian Federation": [55.756, 37.632],
                      "Kosovo": [42.582, 20.902],
                      "Italy": [43.428713, 12.156173],
                      "Greece": [39.243882, 22.078545],
                      "Montenegro": [42.817802, 19.232885],
                      "Croatia": [45.374080, 15.679719],
                      "Serbia": [44.275307, 20.766934],
                      "Austria": [47.427795, 14.240707],
                      "Ukraine": [49.448715, 31.816937],
                      "Slovakia": [48.744892, 19.406103],
                      "Czech Republic": [49.789797, 15.273722],
                      "Ireland": [53.184083, -7.715391],
                      "Denmark": [55.995352, 9.981491],
                      "Slovenia": [46.073618, 14.679797],
                      "Romania": [45.985944, 25.244631],
                      "Albania": [41.147472, 20.072507],
                      "Lithuania": [55.274992, 24.091709],
                      "France": [46.967933, 2.793819],
                      "Sweden": [58.136089, 14.900652],
                      "Turkey": [39.120936, 33.566343],
                      "Switzerland": [46.818993, 8.152376],
                      "Norway": [60.331230, 8.557697],
                      "Hungary": [47.163164, 19.291062],
                      "Netherlands": [52.265760, 5.629735],
                      "Iraq": [33.341973, 44.320839],
                      "India": [23.488056, 78.993478],
                      "Lebanon":[33.890525, 35.830802],
                      "Kuwait":[29.317687, 47.747070],
                      "Saudi Arabia":[24.753444, 46.669013],
                      "Iran":[33.010301, 53.652005],
                      "Yemen":[15.620167, 47.598202],
                      "Bahrain":[26.092226, 50.553897],
                      "Qatar":[25.365382, 51.215336],
                      "Syrian Arab Republic":[35.154318, 38.582985],
                      "Portugal": [39.452106, -8.323991],
                      "United Arab Emirates":[24.022526, 54.737414],
                      "Sri Lanka":[7.524123, 80.842900]}

    countryCenters["Colombia"] = countryCenters["Spain"]                      
    countryCenters["Mexico"] = countryCenters["Spain"]
    countryCenters["Brazil"] = countryCenters["Portugal"]
    countryCenters["Dominican Republic"] = countryCenters["Spain"]
    countryCenters["Norman"] = [49.073779, 0.214844]
    countryCenters["Sardinia"] = [40.023970, 9.075045]
    countryCenters["Tuscany"] = [43.478, 11.149]
        
    passedHeader = False
    for row in kitRows:
        if passedHeader:            
            cell = row[2]
            subGroup = None
            negatives = []
            if cell[0].text is not None:
                subGroup = parseSubgroupText(cell[0].text)
                negatives = parseSubgroupNegatives(cell[0].text)            
                
                kitNo = row[0][0].text          
    
                latCell = row[9][0]
                lonCell = row[10][0]
                lat = None
                lon = None
                if latCell.text is not None:                                
                    lat = float(latCell.text)
                    lon = float(lonCell.text)                
                
                countryName = row[5][0].text
                
                def replaceLocationIfIdInOverrides(latLon, iD):                    
                    for i in idRegionOverrides:
                        if iD == i["id"]:
                            if i["replacement"] in countryCenters:
                                print('replaced', latLon, 'with', countryCenters[i["replacement"]], iD)
                                return countryCenters[i["replacement"]]
                    return latLon
                    
                latLon = [0.0, 0.0]    
                
                if lat is not None and not (lat == 0.0 and lon == 0.0):
                    latLon = [lat, lon]
                else:                                        
                    if countryName in countryCenters:
                        latLon = countryCenters[countryName]                        
                    else:
                        unknownLocations.add(countryName)
                latLon = replaceLocationIfIdInOverrides(latLon, kitNo)
                if latLon[0] != 0.0 and latLon[1] != 0.0:
                    validKits.append({'id': kitNo, 'clade' : subGroup, 'negatives' : negatives, 'latLon': replaceLocationIfIdInOverrides(latLon, kitNo), 'country': countryName})                                        
                
                if subGroup in parseReport:
                    parseReport[subGroup].append(cell[0].text)
                else:
                    parseReport[subGroup] = [cell[0].text]
        else:
            passedHeader = True

    return validKits

def writeFile(kits):
    maskedWriter = open("anonymousKits.txt","w+")
    writer = open("kits.txt","w+")    
    i = 0
    for kit in kits:
        maskedOut = [i, kit["latLon"][0], kit["latLon"][1], kit["clade"], ":".join(kit["negatives"])]
        unmaskedOut = [kit["id"], kit["latLon"][0], kit["latLon"][1], kit["clade"], ":".join(kit["negatives"])]        
        maskedWriter.write(",".join(map(str, maskedOut)))
        maskedWriter.write("\n")
        writer.write(",".join(map(str, unmaskedOut)))
        writer.write("\n")
        i = i + 1

validKits = getValidKits()
writeFile(validKits)
