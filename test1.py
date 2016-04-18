from numpy import *
from pandas import *
from rpy2.robjects.packages import importr
import rpy2.robjects as ro
import pandas.rpy.common as com

ro.r('''                       
newDef <- function(a,b){
   x = runif(10,a,b)
   mean(x)
}   
''')

# calling the newly defined function
r_newDef = ro.globalenv['newDef']
print(r_newDef(4,5))

# calling a builtin function
r_sum = ro.r['sum']
print("sum is", r_sum(ro.IntVector([1,2,3]))[0]) # everything is vector, even what seems like it should be scalar (sum)

