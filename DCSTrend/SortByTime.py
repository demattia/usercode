#!/usr/env/python

def sortContents(contents):
    # Build up a dictionary of timestamps
    contentsDict = {}
    for item in contents:
        print "item =", item
        contentsDict[item] = item.split(",")[0]
        # print contentsDict[item]
    # Sort keys, based on time stamps
    items = contentsDict.keys()
    items.sort(lambda x,y: cmp(contentsDict[x],contentsDict[y]))
    # Report objects in order
    return items


print "Sorting by time"

file = open("full_run.js")



contentsDict = {}

for line in file:
    if line.find("data") != -1:
        array = line.split("[")
        for element in array:
            if len(element) > 0:
                items = element.split("]")[0]
                for item in items.split("\n"):
                    print "time =", item.split(",")[0]
                    contentsDict[item] = item.split(",")[0]
keys = contentsDict.keys()
# print keys
keys.sort(lambda x,y: cmp(contentsDict[x],contentsDict[y]))
print keys
