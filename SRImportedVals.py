from pysr import *
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

Vrv = pd.read_csv("VrvVals.csv")
Vlv = pd.read_csv("VlvVals.csv")
Vspt = pd.read_csv("VsptVals.csv")
tsteps = np.linspace(0,10,num=1000)

xvals = np.column_stack([Vrv,Vlv,tsteps])

SR_model = PySRRegressor(binary_operators = ["+","-","*"],
                         unary_operators = ["cos", "sin", "exp"],
                         elementwise_loss = "f(x, y) = (x-y)^2",
                         model_selection = "score",
                         verbosity = 0)

for i in range(10):
    SR_fit = SR_model.fit(xvals, Vspt.iloc[i,:])
    print(str(SR_fit.sympy()))
