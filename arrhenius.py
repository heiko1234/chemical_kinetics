

import pandas as pd
import numpy as np
from math import e

import plotly.express as px
import plotly.graph_objects as go


from sklearn.linear_model import LinearRegression


# EA = activation Energy (Unit: kJ·mol−1),
# R = universal Gas constante (8,314 J·K−1·mol−1),
# T = absolute (thermodynamic) Temperature (Unit: K, Temp [°C] + 273.15)



# generate data

def k_function(EA, T, A):
    """[summary]

    Args:
        EA ([numeric]): Activation Energy [kJ / mol]
        T ([numeric]): Temperatur [K]
        A ([numeric]): Stossfactor

    Returns:
        [numeric]: activation factor: k [time-1]
    """
    power = -EA*1000/(8.314 * (T+273.15))
    return A * e**power


def calc_data(data, Temp = "Temp", activity="k-value"):
    data["temp_kelvin"]=data[Temp]+273.15
    data["temp_kelvin-1"] = 1/data["temp_kelvin"]
    data["ln_k"] = np.log(abs(data[activity]))
    return data


def plot_kinetics(data, Temp = "Temp", activity="k-value", color = None):

    data = calc_data(data=data, Temp = Temp, activity=activity)

    fig = px.scatter(data, x = "Temp-1", y = "ln_k", color=color)
    fig.show()


def get_EA(data, Temp_1 = "temp_kelvin-1", ln_k="ln_k"):

    X = np.array(data[Temp_1])
    X = X.reshape(-1,1)
    y = np.array(data[ln_k])

    reg = LinearRegression().fit(X, y)
    # reg.score(X, y)

    return (round(abs(reg.coef_[0]*8.314)/1000, 3) )


def get_A(data, EA, Temp = "Temp", kvalue="k-value"):
    # EA = 45
    # Temp = "Temp"
    # kvalue="k-value"

    power = -EA*1000/(8.314 * (data[Temp]+273.15))
    exp_value = e**power
    exp_value

    X = np.array(exp_value)
    X = X.reshape(-1,1)
    y = np.array(data[kvalue])

    reg = LinearRegression().fit(X, y)
    reg.coef_
    reg.intercept_

    return (abs(reg.coef_[0]))


def translate(data, columnlist):
    for element in columnlist:
        try:
            data[element] = data[element].str.replace(",", ".")
            data[element] = pd.to_numeric(data[element], errors='coerce')
        except:
            continue
    return data



def kinetic_1_order(c0, k, t):
    c = c0 * e**(-k*t)
    return c



def group_analysis(data, column = "name", filter="Ord_1", vary="temp"):
    #column= "name"
    #filter = "Ord_1"
    #vary= "temp"

    idata = data[data[column] == filter]
    idata = idata.reset_index(drop = True)
    idata

    timestamps = idata.iloc[:,4:].columns

    timestamps = list(timestamps)
    timestamps = [int(i) for i in timestamps]

    timestamps

    output = pd.DataFrame(columns = [])

    output_vary=[]
    output_column=[]
    output_value=[]

    for i in range(idata.shape[0]):

        cname = idata[vary][i]

        jdata = idata.iloc[i:i+1,4:]

        analysis_values = jdata.iloc[0,:]
        analysis_values = list(analysis_values)

        analysis_values
        analysis_values_ln = [np.log(i/analysis_values[0]) for i in analysis_values]
        analysis_values_ln

        df = pd.DataFrame()
        df["x"] = timestamps
        df["y"] = analysis_values_ln
        df
        df_modelling = df.dropna()
        df_modelling = df_modelling.reset_index(drop = True)


        # regression
        reg = LinearRegression().fit(np.vstack(df_modelling['x']), df_modelling["y"])
        df['bestfit'] = reg.predict(np.vstack(df['x']))
        df

        # model
        reg.coef_[0]

        # output[cname] = [reg.coef_[0]]

        output_column.append(filter)
        output_vary.append(cname)
        output_value.append(abs(reg.coef_[0]))
    
    output[column] = output_column
    output[vary] = output_vary
    output["value"] = output_value
    
    return output


def table_analysis(data, column= "name", names=["Ord_1"], vary = "temp"):
    output = []
    for i in names:
        df = group_analysis(data=data, column = column, filter=i, vary=vary)
        output.append(df)
    
    df_output = pd.concat(output, axis =0)
    df_output = df_output.reset_index(drop = True)

    return df_output




# kinetic 1. order

# c[A] = c0 * exp(-k*t)

# t = time
# k = activity coefficient
# c0 = Conentration at the beginning, at time t = 0

# ln (c[A] / c0) = -k * t

data = pd.read_csv("/home/heiko/Repos/chemical_kinetics/sample_csv.csv")

data.columns
list(data.columns).index("0") #4

data = translate(data, columnlist=['k', '0', '10', '20', '30', '40', '60', '80','100', '160'])


data

data[data["name"] == "Ord_1"]
data[data["temp"] == 20]


list(set(data["name"]))



group_analysis(data = data, column = "name", filter="Ord_1", vary="temp")

group_analysis(data = data, column = "name", filter="Exp2", vary="temp")

group_analysis(data = data, column = "temp", filter=20, vary="name")


gdata = group_analysis(data = data, column = "name", filter="Ord_1", vary="temp")
gdata

cdata = calc_data(data=gdata, Temp = "temp", activity="value")
cdata


tdata = table_analysis(data, column= "name", names=["Ord_1"], vary = "temp")
tdata


list_names= list(set(data["name"]))  # all unique elements in name get calculated

tdata = table_analysis(data, column= "name", names=list_names, vary = "temp")
tdata

cdata = calc_data(data=tdata, Temp = "temp", activity="value")
cdata


get_EA(cdata, Temp_1 = "temp_kelvin-1", ln_k="ln_k")   #50.948

get_A(cdata, EA = 50.948, Temp = "temp", kvalue="value")  #11,972,690


k_function(EA= 50.948, T = 30, A= 11972690)

kinetic_1_order(c0 = 100, k = 0.009999, t = 10)
kinetic_1_order(c0 = 90.4, k = 0.02, t = 10)


fig = px.scatter(data_frame = cdata, x = "temp_kelvin-1", y="ln_k", color="name")
fig.show()





###

timestamps = data.iloc[0:1,4:].columns
timestamps = list(timestamps)
timestamps = [int(i) for i in timestamps]

timestamps

i = 3

analysis_values = data.iloc[i:i+1, 4:]
analysis_values = analysis_values.iloc[0,:]
analysis_values = list(analysis_values)

analysis_values
analysis_values_ln = [np.log(i/analysis_values[0]) for i in analysis_values]
analysis_values_ln


fig = px.scatter(x = timestamps, y = analysis_values)
fig.show()

# fig = px.scatter(x=timestamps, y=analysis_values_ln)
# fig.show()

df = pd.DataFrame()
df["x"] = timestamps
df["y"] = analysis_values_ln
df
df_modelling = df.dropna()
df_modelling = df_modelling.reset_index(drop = True)

# regression
reg = LinearRegression().fit(np.vstack(df_modelling['x']), df_modelling["y"])
df['bestfit'] = reg.predict(np.vstack(df['x']))
df

# plotly figure setup
fig=go.Figure()
fig.add_trace(go.Scatter(name='X vs Y', x=df['x'], y=df['y'].values, mode='markers'))
fig.add_trace(go.Scatter(name='line of best fit', x=df["x"], y=df['bestfit'], mode='lines'))

# plotly figure layout
fig.update_layout(xaxis_title = 'X', yaxis_title = 'Y')

fig.show()


# model
reg.coef_[0]
reg.intercept_
reg.score(np.vstack(df_modelling['x']), df_modelling["y"]) 


# generate some k values

k_function(EA= 45, T = 50, A= 40000)  #39335.6
k_function(EA= 45, T = 40, A= 40000)  #39314.6
k_function(EA= 45, T = 30, A= 40000)  #39292.2
k_function(EA= 45, T = 20, A= 40000)  #39268.2


# make dataframe

data = pd.DataFrame()
data["Temp"] = [20, 30, 40, 50]
data["k-value"] = [39268, 39292, 39315, 39335]
data["name"] = "Test1"

data2 = pd.DataFrame()
data2["Temp"] = [20, 30, 40, 50]
data2["k-value"] = [39266, 39294, 39314, 39333]
data2["name"] = "Test2"

data3 = pd.DataFrame()
data3["Temp"] = [21, 31, 41, 49]
data3["k-value"] = [39270, 39300, 39310, 39328]
data3["name"] = "Test3"



# concat the data

data = pd.concat([data,data2,data3], axis=0)
data


# analyse the data

data = calc_data(data=data, Temp = "Temp", activity="k-value")
data


plot_kinetics(data, color="name")

data

get_EA(data, Temp_1 = "Temp-1", ln_k="ln_k")

get_A(data, EA = get_EA(data), Temp = "Temp", kvalue="k-value")


# test these values

data_test = pd.DataFrame()
data_test["Temp"] = [20, 30, 40, 50]
data_test["k-value"] = k_function(EA= 45, T=data_test["Temp"], A= 40050)
data_test
data_test["name"] = "testdata"



plot_kinetics(data=data_test, Temp="Temp", activity="k_values")



data

data = pd.concat([data, data_test], axis=0)

data


plot_kinetics(data, Temp = "Temp", activity="k-value", color = "name")





