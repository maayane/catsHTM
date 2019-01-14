import numpy as np

ID=[[[[[[3833.0]]], [[[3913.0]]], [[[3961.0]]]], [[[[4553.0]]], [[[4665.0]]], [[[4745.0]]]]]]

print('the nested list is',ID)

def simplify_list(val):
        if isinstance(val, list) == False:
            return val
        else:
            if len(val) > 1:
                return val
            else:
                return simplify_list(val[0])

'''IDx=[]
for i in ID:
    a=simplify_list(i)
    if isinstance(a, (list, tuple, np.ndarray))==False:
        IDx.append(simplify_list(i))
    else:
       for j in a:
           IDx.append(simplify_list(j))

print('IDx is',IDx)
'''
def simplify2(x):
    IDc=[]
    for i in x:
        if isinstance(i, (list, tuple, np.ndarray)) == True:
            for j in i:
                IDc.append(j)
        else:
            IDc.append(i)
    return IDc
    #return simplify2(IDc)

def simplify3(x):
    if isinstance(x[0],(list, tuple, np.ndarray)) == False:
        return x
    else:
        y=simplify2(x)
        print(y)
        return simplify3(y)



print ('simplify3(IDx) is',simplify3(ID))

'''

IDy=[]
for i in IDx:
    if isinstance(i, (list, tuple, np.ndarray)) == True:
        for j in i:
            IDy.append(j)
    else:
        IDy.append(i)
IDz=[]
for i in IDy:
    if isinstance(i, (list, tuple, np.ndarray)) == True:
        for j in i:
            IDz.append(j)
    else:
        IDz.append(i)

IDw=[]
for i in IDz:
    if isinstance(i, (list, tuple, np.ndarray)) == True:
        for j in i:
            IDw.append(j)
    else:
        IDw.append(i)
IDh=[]
for i in IDw:
    if isinstance(i, (list, tuple, np.ndarray)) == True:
        for j in i:
            IDh.append(j)
    else:
        IDh.append(i)
'''

print('the list without brakets is',simplify3(ID))
