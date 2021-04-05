
In [40]: a=np.array(a)                                                                                      

In [41]: a                                                                                                  
Out[41]: 
array([[0, 2],
       [0, 3],
       [2, 3]])

In [42]: b                                                                                                  
Out[42]: 
array([[False, False, False, False],
       [False, False, False, False],
       [False, False, False, False],
       [False, False, False, False]])

In [43]: c0=a[:,0]                                                                                          

In [44]: c1=a[:,1]                                                                                          

In [45]: b[c0,c1] = 1                                                                                       

In [46]: b                                                                                                  
Out[46]: 
array([[False, False,  True,  True],
       [False, False, False, False],
       [False, False, False,  True],
       [False, False, False, False]])
