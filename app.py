#!/usr/bin/env python

#importing required libraries


from io import StringIO 
import pandas as pd
import numpy as np
import tensorflow as tf
import time
from collections import Counter
import matplotlib.pyplot as plt
import re
import io
import plotly.graph_objects as go
import base64
import os
import joblib
import csv
from clocks import han2020_coefs
from sklearn.linear_model import LinearRegression
from sklearn.impute import SimpleImputer
from sklearn.experimental import enable_iterative_imputer
from sklearn.impute import IterativeImputer
import math


AltumAge_cpgs = np.array(pd.read_pickle('clocks/altum_cpgs.pkl'))
scaler = pd.read_pickle('clocks/altum_scaler.pkl')
AltumAge = tf.keras.models.load_model('clocks/AltumAge.h5')

#all_grim_cpgs = grim_cpgs.iloc[2:, 1].tolist()

df_cpg_coeff = pd.read_csv("clocks/hannum_coeffs.csv").set_index("Marker")

#all_hannum_cpgs = df_cpg_coeff.index.tolist()

df_slope = pd.read_csv("clocks/horvath_coeffs.csv")

#all_horvath_cpgs = df_slope.iloc[1:, 0].tolist()

df_coeff_pheno = pd.read_csv("clocks/pheno_coeffs.csv").set_index("CpG")
df_horvath_sb = pd.read_csv("clocks/horvath_skin_blood_coeff.csv")


zhangen_params = pd.read_table('clocks/zhangen.coef', sep=' ')
zhangen_params['feature'] = zhangen_params['probe']
zhangen_params['coefficient'] = zhangen_params['coef']
zhangen_cpgs = zhangen_params['feature'][1:].tolist()

yingCausAge_params = pd.read_csv('clocks/YingCausAge.csv')
yingCausAge_params['feature'] = yingCausAge_params['term']
yingCausAge_params['coefficient'] = yingCausAge_params['estimate']
yingCausAge_cpgs = yingCausAge_params['feature'][1:].tolist()

han2020_cpgs = han2020_coefs.cpgs


import streamlit as st
st.title("Ensemble Aging Clock :alarm_clock:")


input_file = st.file_uploader("Please choose a CSV file; Download example file here: [example](https://github.com/hayanlee/EnsembleMeAgingClock/blob/main/examples/individual01.csv)")




# Check if a file has been uploaded
#if input_file is not None:
    # Process the uploaded file and display a message
 #   st.write("File uploaded and processed.")'''

#input_file = st.file_uploader("Please choose a CSV file [Single Sample Example File](examples/individual01.csv) [Multi Sample Example File](examples/dat0BloodIllumina450K.csv)")



#selected_tab = st.selectbox("Select Sample Type", ["Single Sample", "Multiple Sample", "Single Sample With GrimAge"])

def get_altum_age(df_meth):
    df_meth = df_meth.set_index(df_meth.columns[0])
    meth_data = df_meth.loc[AltumAge_cpgs].transpose()
    meth_scaled = scaler.transform(meth_data)
    pred_age_AltumAge = AltumAge.predict(meth_scaled).flatten()
    return pred_age_AltumAge

def get_hannum_age(df_meth):
    #print('hannum1', file=sys.stderr)
    df_meth = df_meth.set_index(df_meth.columns[0])
    
    #print('hannum2', file=sys.stderr)
    # Get the common row labels between df1 and df2
    common_labels = df_cpg_coeff.index.intersection(df_meth.index)
    # Select only the rows in df2 that have matching row labels with df1 in the same order
    df_meth = df_meth.loc[common_labels]
    #print('hannum3', file=sys.stderr)
    coef = np.array(df_cpg_coeff["Coefficient"][1:] )
    #coef
     
    intercept = 18.695467523254287
    
    
    sample_methylation = np.array(df_meth[1:])
    #print('hannum', file=sys.stderr)
    #print(sample_methylation, file=sys.stderr)
    #print(sample_methylation.shape, file=sys.stderr)
    #print(len(coef), file=sys.stderr)
    pred = np.matmul(coef, sample_methylation) - intercept
    #print(pred, file=sys.stderr)
    #print('hannum4', file=sys.stderr)
    return pred



def get_horvath_age(df_cpg):
    df_cpg = df_cpg.set_index(df_cpg.columns[0])

    df_slope_only = df_slope.iloc[0:]
    df_slope_only = df_slope_only.iloc[:, :2]
    df_slope_only = df_slope_only.set_index("CpGmarker")

    coef = np.array( df_slope_only['CoefficientTraining'][1:] )


    input2 = df_cpg.loc[df_slope_only.index[1:], :].T
    
    pred = np.matmul (input2, coef) + df_slope_only['CoefficientTraining'][0]

    def anti_transform_age(exps):
        import numpy as np
        adult_age = 20
        ages = []
        for exp in exps:
            import numpy as np
            if exp < 0:
                age = (1 + adult_age)*(np.exp(exp))-1
                ages.append(age)
            else:
                age = (1 + adult_age)*exp + adult_age
                ages.append(age)
        ages = np.array(ages)
        return ages

    age = anti_transform_age(pred)
    return age


def get_pheno_age(df):
    
    global df_coeff_pheno
    df_cpg_coeff_pheno = df_coeff_pheno.iloc[:, [4]]
    intercept = df_cpg_coeff_pheno.iloc[0,0]
    df_cpg_coeff_pheno = df_cpg_coeff_pheno.drop(df_cpg_coeff_pheno.index[0])

    
    df_sample = df.set_index(df.columns[0])
    
   
    # Get the common row labels between df1 and df2
    common_labels = df_cpg_coeff_pheno.index.intersection(df_sample.index)
    # Select only the rows in df2 that have matching row labels with df1 in the same order
    df_sample = df_sample.loc[common_labels]
    
    coef = np.array(df_cpg_coeff_pheno["Weight"])
    
    
    sample_methylation = np.array(df_sample)
    


    pred = np.matmul(coef, sample_methylation) + intercept
  
    return pred

def get_horvath_sb_age(df):
    
    global df_horvath_sb
    
    df_cpg_coeff_horvath_sb = df_horvath_sb.set_index(df_horvath_sb.columns[0])
    df_cpg_coeff_horvath_sb = df_cpg_coeff_horvath_sb.iloc[:, 0]

    intercept = df_cpg_coeff_horvath_sb.iloc[0]
    

    df_sample = df.set_index(df.columns[0])
    
    #print(df_sample)
   

    # Get the common row labels between df1 and df2
    common_labels = df_cpg_coeff_horvath_sb.index.intersection(df_sample.index)
    
    # Select only the rows in df2 that have matching row labels with df1 in the same order
    df_sample = df_sample.loc[common_labels]
    
    

    coef = np.array(df_cpg_coeff_horvath_sb[1:])
    coef = np.ravel(coef)


    sample_methylation = np.array(df_sample)
    #sample_methylation = np.ravel(sample_methylation)

    
    
    pred = np.matmul(sample_methylation.T, coef) + intercept

    def anti_transform_age(exps):
            import numpy as np
            adult_age = 20
            ages = []
            for exp in exps:
                import numpy as np
                if exp < 0:
                    age = (1 + adult_age)*(np.exp(exp))-1
                    ages.append(age)
                else:
                    age = (1 + adult_age)*exp + adult_age
                    ages.append(age)
            ages = np.array(ages)
            return ages

    age = anti_transform_age(pred)
    return age 

def get_ZhangEn_age(df):
    def scale_row(x):
        """
        Scales the input data per row with mean 0 and std 1.
        """
        row_means = np.mean(x, axis=1, keepdims=True)
        row_stds = np.std(x, axis=1, keepdims=True)

        # Avoid division by zero in case of a row with constant value
        row_stds[row_stds == 0] = 1

        x_scaled = (x - row_means) / row_stds
        return x_scaled
    df_cpg = df.set_index(df.columns[0])
    weights = np.array(zhangen_params['coefficient'][1:].tolist())
    intercept = np.array([zhangen_params['coefficient'][0]])
    X = df_cpg.loc[zhangen_cpgs].values.T
    X = scale_row(X)
    Y = np.dot(X, weights) + intercept
    return Y
    
def get_Han2020_age(df):
    df_cpg = df.set_index(df.columns[0])
    X_data = df_cpg.loc[han2020_coefs.cpgs].values.T
    weights = np.array(han2020_coefs.coefficients[1:])
    intercept = np.array(han2020_coefs.coefficients[0])
    Y = np.dot(X_data, weights) + intercept
    def anti_log_linear(x, adult_age=20):
        """
        Applies an anti-logarithmic linear transformation to a value.
        """

        exp_transform = (1 + adult_age) * np.exp(x) - 1
        lin_transform =  (1 + adult_age) * x + adult_age
        result = np.where(x < 0, exp_transform, lin_transform)
        return result
    return anti_log_linear(Y)

def get_YingCaus_age(df):
    df_cpg = df.set_index(df.columns[0])
    weights = np.array(yingCausAge_params['coefficient'][1:].tolist())
    intercept = np.array([yingCausAge_params['coefficient'][0]])
    X = df_cpg.loc[yingCausAge_cpgs].values.T
    Y = np.dot(X, weights) + intercept
    return Y

def is_number(string):
    try:
        float(string)  # Attempt to convert the string to a float
        return True
    except ValueError:
        return False

st.write("")
st.write("Select clocks for EnsembleNaive. The default selection has been shown to be the most optimal.")


st.write("")

# Create columns to display checkboxes side by side
col1, col2, col3, col4 = st.columns(4)

# Checkbox for each clock model in respective columns

with col1:
    
    altum_selected = st.checkbox('Altum', value=True)
    hannum_selected = st.checkbox('Hannum', value=False)
    

with col2:
    horvath_selected = st.checkbox('Horvath', value=True)
    pheno_selected = st.checkbox('Pheno', value=False)
    

with col3:
    horvath_sb_selected = st.checkbox('Skin Blood', value=True)
    zhangEn_selected = st.checkbox('Zhang2019', value=True)
    
with col4:
    han2020_selected = st.checkbox('Han2020', value=True)
    yingCausal_selected = st.checkbox('Ying Causal', value=True)


col1, col2 = st.columns([3, 1])


if input_file is not None:
    print("1")
    # Function to find missing values and compare with cpg lists
    
    
    df = pd.read_csv(input_file, sep=r'[,\t]', engine='python')
    file_list = ['tools/altum_cpgs_list.txt', 'tools/hannum_cpgs_list.txt', 'tools/horvath_cpgs_list.txt', 'tools/pheno_cpgs_list.txt', 'tools/horvath_SB_cpgs_list.txt', 'tools/han_cpgs_list.txt', 'tools/ying_cpgs_list.txt', 'tools/zhang_cpgs_list.txt']
    # Create an empty list to store the overlapping values
    overlapping_cpgs = []
    all_missing_values = []
    # Iterate through the files in file_list
    for file in file_list:
        count = 0
        with open(file, 'r') as file:
            cpgs_list = file.readlines()

        cpgs_list = [line.strip() for line in cpgs_list]
    

        # Filter the DataFrame rows to only include rows with overlapping values
        overlapping_values = df[df.iloc[:, 0].isin(cpgs_list)]

        # Extend the overlapping_cpgs list with the filtered values from the DataFrame
        overlapping_cpgs.extend(overlapping_values.values.tolist())

        # Convert the overlapping_cpgs list into a DataFrame
        overlapping_df = pd.DataFrame(overlapping_cpgs, columns=df.columns)

        # Now you can iterate over this overlapping DataFrame as you did before
        for i, row in overlapping_df.iloc[1:].iterrows():
            for j, val in enumerate(row[1:], start=1):
                if pd.isnull(val) or val == 0:
                    count = count + 1
                    
                    #print("For" , file, ": ", "Row Label:", overlapping_df.iloc[i, 0], "Column Label:", overlapping_df.columns[j])
        
        file_name = file.name
        print("For" , file_name.replace("tools/", "").replace("_cpgs_list.txt", ""), ": there are: ", count, " 0 or NA values")
        
        missing_values = set(cpgs_list) - set(df.iloc[:, 0])
        all_missing_values = all_missing_values + list(missing_values)
 
        print("For" , file_name.replace("tools/", "").replace("_cpgs_list.txt", ""), ": there are: ", len(missing_values), "missing CpGs")
    
    all_missing_values = list(set(all_missing_values))
    for value in all_missing_values:
        new_row = [value] + [np.nan] * (df.shape[1] - 1)  # Create a new row with value and 0.5 in other columns
        # print(new_row)
        df.loc[len(df)] = new_row  # Add new row to DataFrame
    
    labels = df.iloc[:, 0]  # Assuming the first column contains the CPG labels
    methylation_data = df.iloc[:, 1:]  # Assuming the methylation data starts from the second column

    df.fillna(np.nan, inplace=True)

    # Impute missing values using IterativeImputer
    imputer_iterative = IterativeImputer(max_iter=10, random_state=0)  # You can adjust max_iter as needed
    imputed_data_iterative = imputer_iterative.fit_transform(methylation_data)

    # Convert imputed data back to DataFrame
    # Create the DataFrame
    imputed_df_iterative = pd.DataFrame(imputed_data_iterative, columns=methylation_data.columns)

    # Add the labels as the first column
    imputed_df_iterative.insert(0, 'Labels', labels)
 
    df = imputed_df_iterative


    

    if True:
    
        
        

        # Calculate predictions
        preds = []
    
        altum_prediction = str(round(np.float32(get_altum_age(df)[0]), 2))
        hannum_prediction = str(round(np.float32(get_hannum_age(df)[0]), 2))
        horvath_prediction = str(round(np.float32(get_horvath_age(df)[0]), 2))
        pheno_prediction = str(round(np.float32(get_pheno_age(df)[0]), 2))
        horvath_sb_prediction = str(round(np.float32(get_horvath_sb_age(df)[0]), 2))
        zhangEn_pred = str(round(np.float32(get_ZhangEn_age(df)[0]), 2))
        han2020_pred = str(round(np.float32(get_Han2020_age(df)[0]), 2))
        yingCausal_pred = str(round(np.float32(get_YingCaus_age(df)[0]), 2))
        if altum_selected:
            preds.append(altum_prediction)
        if hannum_selected:
            preds.append(hannum_prediction)
        if horvath_selected:
            preds.append(horvath_prediction)
        if pheno_selected:
            preds.append(pheno_prediction)
        if horvath_sb_selected:
            preds.append(horvath_sb_prediction)
        if zhangEn_selected:
            preds.append(zhangEn_pred)
        if han2020_selected:
            preds.append(han2020_pred)
        if yingCausal_selected:
            preds.append(yingCausal_pred)

        # Array of data
        data = [float(val) for val in preds]
        ensembleLR_pred = round(19.275324793823557 + (0.3437195049507392 * float(altum_prediction)) + (-0.20633851531241806 * float(han2020_pred)) + (-0.5173259165634722 * float(hannum_prediction)) + (0.2836010547077375 * float(horvath_prediction)) + (0.1857810346592119 * float(horvath_sb_prediction)) + (-0.02208310995460718 * float(pheno_prediction)) + (-0.08526975203798473 * float(yingCausal_pred)) + (0.5640149630320757 * float(zhangEn_pred)), 2)
        average = np.average(data)
        average = round(average, 2)
        average = str(np.float32(average))

        preds.append(average)

        with col2: 
            st.write(" ")
            st.write("EnsembleNaive: " + average + "*")
            st.write("EnsembleLR: " + str(ensembleLR_pred) + "*")
            text= []
            if altum_selected:
                st.write("AltumAge: " + altum_prediction)
                text.append("Altum")
            if hannum_selected:
                st.write("Hannum: " + hannum_prediction)
                text.append("Hannum")
            if horvath_selected:
                st.write("Horvath: " + horvath_prediction)
                text.append("Horvath")
            if pheno_selected:
                st.write("PhenoAge: " + pheno_prediction)
                text.append("Pheno")
            if horvath_sb_selected:
                st.write("Skin Blood: " + horvath_sb_prediction)
                text.append("Skin Blood")
            if zhangEn_selected:
                st.write("Zhang2019 Clock: " + zhangEn_pred)
                text.append("Zhang2019")
            if han2020_selected:
                st.write("Han2020 Clock: " + han2020_pred)
                text.append("Han2020")
            if yingCausal_selected:
                st.write("YingCausalAge: " + yingCausal_pred)
                text.append("YingCausal")

            std_dev = np.std(data)
            std_dev = round(std_dev, 2)
            std_dev = str(std_dev)
            st.write("---")
            st.write("Standard Deviation: " + std_dev)


        with col1: 
            fig = go.Figure()
            
            fig.add_trace(go.Scatter(
                x = [-1.0], 
                y = np.array([preds[-1]]),
                mode = 'markers + text',
                marker=dict(size=8, color = 'rgb(255, 0, 0)'),  # Adjust marker size here
                text=['EnsembleNaive: ' + str(preds[-1])],  # Custom hover text
                textposition="middle right",
                textfont=dict(size=10),
                showlegend=False
            ))
            fig.add_trace(go.Scatter(
                x = [-0.9], 
                y = np.array([ensembleLR_pred]),
                mode = 'markers + text',
                marker=dict(size=8, color = 'rgb(255, 0, 0)'),  # Adjust marker size here
                text=['EnsembleLR: ' + str(ensembleLR_pred)],  # Custom hover text
                textposition="middle right",
                textfont=dict(size=10),
                showlegend=False
            ))
            coords = [-0.80, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1]
            custom_colors = ['#AFE1AF', '#00A7E1', '#89CFF0', '#50C878', '#009245', '#F08080', '#FFD700', '#9370DB']
            
            selected_preds = [altum_selected, hannum_selected, horvath_selected, pheno_selected, horvath_sb_selected, zhangEn_selected, han2020_selected, yingCausal_selected]
            
            for i, (pred, selected) in enumerate(zip(preds[:-1], text)):
                    fig.add_trace(go.Scatter(
                        x=[coords[i]],  # x-coordinate for the point
                        y=[pred],  # y-coordinate for the point
                        mode='markers + text',
                        text=[text[i] + ': ' + str(pred)],  # Custom hover text
                        textposition="middle right",
                        textfont=dict(size=10),
                        marker=dict(size=8, color=custom_colors[i]),
                        name=f'Point {i+1}',  # Name for the legend
                        showlegend=False
                    ))
            
            fig.add_trace(go.Box(
                y=preds[:-1],
                name="Clock Predictions",
                jitter=0.3,
                pointpos=-1.8,
                line_color='rgb(161, 236, 254)',
                showlegend=False
            ))

            numeric_preds = [float(value) for value in preds[:-1]]
            max_value = max(numeric_preds) + 5
            min_value = min(numeric_preds) - 5
            fig.update_yaxes(range=[min_value, max_value])
            fig.update_yaxes(tickfont=dict(size=20))
            fig.update_layout(
                title_text="Ensemble Aging Clock Predictions",                          
                shapes=[
                    dict(
                        type="rect",
                        xref="paper",
                        yref="paper",
                        x0=0,
                        y0=0,
                        x1=1,
                        y1=1,
                        line=dict(
                            color="black",
                            width=1
                        ),
                        fillcolor="rgba(0,0,0,0)",
                    )
                ],
                xaxis=dict(
                    showticklabels=False  # Remove x-axis tick labels
                ),
                width=650,  # Adjust width to make the plot narrower
                height=700,  # Adjust height to make the plot longer
            )

            st.plotly_chart(fig, use_container_width=True)
                
    

        
            