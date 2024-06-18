import pandas as pd 
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from scipy.stats import t
from scipy.stats import normaltest 
from scipy.stats import chi2_contingency
from scipy.stats import f_oneway


dataframe= "C:/Users/gonza/OneDrive/Escritorio/Doctorado/Genomas/S_meliloti/Sinorhizobium_meliloti.csv"

df= pd.read_csv (dataframe)

def functions (df):
    top_categories= df["product"].value_counts().head(10)
    other_categories= df["product"].value_counts().tail(-10).sum()
    df_combined_categories= pd.concat ([top_categories, pd.Series ({"Others":other_categories})])
    colores = ['lightcoral', 'lightskyblue', 'lightgreen', 'lightpink', 'lightsalmon', 'lightyellow', 'lightblue', 'lightcyan', 'lightgray', 'orange', 'lightseagreen']
    pie, _ = plt.pie(df_combined_categories, startangle=90, colors=colores)

    porcentajes = (df_combined_categories / df_combined_categories.sum() * 100)
    labels = [f'{categoria} - {porcentaje:.1f}%' for categoria, porcentaje in zip(df_combined_categories.index, porcentajes)]

    plt.legend(pie, labels, loc="upper left", bbox_to_anchor= (1,1))

    plt.savefig ("products_smorfs_Sme", bbox_inches= "tight")
    plt.close ()

def histogram (df):
    len_column = df["len"]
    in_order= len_column.sort_values ()
    range_bins = [10, 20, 30, 40, 50, 60, 70]

    sns.displot (len_column, bins=range_bins, kde=False)
    plt.xticks ()
    plt.ylabel ("Frecuencies")
    plt.xlabel ("smORFS size")
    plt.savefig ("histogram_smORFs_Sme")
    plt.close ()

def boxplot(df):
    unique_genomes = df['genome'].unique()

    bins_list = []
    count_list = []

    for genome in unique_genomes:
        subset_df = df[df['genome'] == genome]

        subset_df["bins"] = pd.cut(subset_df["len"], bins=[0, 10, 20, 30, 40, 50, 60, 70])

        bins = subset_df.groupby("bins")["len"].count().reset_index()
        bins_list.extend(bins['bins'])
        count_list.extend(bins['len'])

    new_df = pd.DataFrame({'bins': bins_list, 'count': count_list})

    sns.boxplot(data=new_df, x="bins", y="count")
    plt.xlabel("SEPs size (aa)")
    plt.ylabel("Count")
    plt.xticks(rotation=45)
    plt.savefig("boxplot_smorfs_Sme", bbox_inches= "tight")
    plt.close()

def piechart (df):

    priority_order = ["Conserved", "Core", "Soft core", "Shell", "Cloud", "Probably not conserved"]
    df['Results_column'] = pd.Categorical(df['Results_column'], categories=priority_order, ordered=True)

    classification= df.groupby ("Results_column") ["genome"].count ().reset_index ()
    frecuencies= classification.iloc [:,1]
    category= classification.iloc [:, 0]
    colors = ["deepskyblue", "lightblue", "lightgreen", "violet", "lightcoral", "crimson"]

    pie, _ = plt.pie(frecuencies, startangle=90, colors=colors)

    porcent= (classification["genome"]/classification ["genome"].sum ()*100)
    labels = ['{0} - {1:1.2f} %'.format(i,j) for i,j in zip(category, porcent)]

    plt.legend(pie, labels, loc="upper left", bbox_to_anchor= (1,1))
    plt.savefig ("smORFs_classified_Sme.png", bbox_inches= "tight")


functions (df)
histogram (df)
boxplot (df)
piechart (df)

#   MEDIDAS DE CENTRALIZACION
media = df["len"].mean(axis=0, skipna=True, numeric_only=False)

print ("La media es " + str (media))

mediana = df["len"].median(axis=0, skipna=True, numeric_only=False)

print ("La mediana es " + str (mediana))

moda = stats.mode(df["len"], nan_policy='omit')
print (f"la moda es: {moda}")

k= 0.5
percentil_k = df["len"].quantile(q=k)
print("El percentil " + str(k*100) + " es " + str(percentil_k))

#   MEDIDAS DE DISPERSION

rango = df["len"].max()-df["len"].min()
RI = df["len"].quantile(0.75)-df["len"].quantile(0.25)

print ("El rango es " + str (rango))
print ("El recorrido intercuartílico es " + str (RI))

desviacion_tipica = df["len"].std(axis=None, skipna=True, ddof=1, numeric_only=False)
varianza = df["len"].var(axis=None, skipna=True, ddof=1, numeric_only=False)

print ("La desviación típica es " + str (desviacion_tipica))
print ("La varianza es " + str (varianza))

#   ASIMETRIA

#Fisher

skewness = df["len"].skew(axis=0, skipna=True, numeric_only=False)
print (skewness)

# Curtosis

kurtosis = df["len"].kurt(axis=0, skipna=True, numeric_only=False)
print (kurtosis)


#   Promedio de smORFs por genoma.

df_dict = {}

lengths = []

for genoma in df['genome'].unique():
    df_dict[genoma] = df[df['genome'] == genoma]
    
    
    lengths.append(len(df_dict[genoma]))

mean_length = pd.Series(lengths).mean()
print(f"La media del promedio de smORFS por genoma : {mean_length}")

mediana_lenght= pd.Series(lengths).median()
print (f"La mediana del promedio de smORFs por genoma es: {mediana_lenght}")

desviacion_tipica = pd.Series(lengths).std(axis=None, skipna=True, ddof=1, numeric_only=False)
varianza = pd.Series(lengths).var(ddof=1)

print ("La desviación típica del promedio de smORFS por genoma es: " + str (desviacion_tipica))
print ("La varianza del promedio de smORFs por genoma es: " + str (varianza))

tamaño_muestra = len (df_dict)

nivel_confianza = 0.95

grados_libertad = tamaño_muestra - 1

valor_critico = t.ppf((1 + nivel_confianza) / 2, df=grados_libertad)

error_estandar_media = desviacion_tipica / (tamaño_muestra ** 0.5)

#límites del intervalo
limite_inferior = mean_length - (valor_critico * error_estandar_media)
limite_superior = mean_length + (valor_critico * error_estandar_media)

print(f"Intervalo de confianza al {nivel_confianza * 100}%: ({limite_inferior}, {limite_superior})")

#   Genes Conservados
priority_order = ["Conserved", "Core", "Soft core", "Shell", "Cloud", "Probably not conserved"]

df["Results_column"] = pd.Categorical(df["Results_column"], categories=priority_order, ordered=True)

genes50 = df[df["Results_column"] < "Shell"]

df_dict_conservados = {}

lengths_conservados = []

for genoma in genes50['genome'].unique():
    df_dict[genoma] = genes50[genes50['genome'] == genoma]
    
    
    lengths_conservados.append(len(df_dict[genoma]))


mean_length_conservados = pd.Series(lengths_conservados).mean()
print(f"La media de smORFs conservados por genoma es: {mean_length_conservados}")

mediana_lenght_conservados= pd.Series(lengths_conservados).median()
print (f"La mediana del promedio de smORFs conservados por genoma es: {mediana_lenght_conservados}")

desviacion_tipica_conservados = pd.Series(lengths_conservados).std(axis=None, skipna=True, ddof=1, numeric_only=False)
varianza_conservados = pd.Series(lengths_conservados).var(ddof=1)

print ("La desviación típica del promedio de smORFS conservados por genoma es: " + str (desviacion_tipica_conservados))
print ("La varianza del promedio de smORFs conservados por genoma es: " + str (varianza_conservados))

tamaño_muestra = len (df_dict)

nivel_confianza = 0.95

grados_libertad = tamaño_muestra - 1

valor_critico_conservados = t.ppf((1 + nivel_confianza) / 2, df=grados_libertad)

error_estandar_media_conservados = desviacion_tipica_conservados / (tamaño_muestra ** 0.5)

# intervalo de confianza
limite_inferior_conservados = mean_length_conservados - (valor_critico_conservados * error_estandar_media_conservados)
limite_superior_conservados = mean_length_conservados + (valor_critico_conservados * error_estandar_media_conservados)

print(f"Intervalo de confianza al {nivel_confianza * 100}%: ({limite_inferior_conservados}, {limite_superior_conservados}) para smORFs conservados")

normaltest_result = normaltest(df["len"], nan_policy='omit')
print("Normaltest p-value:", normaltest_result.pvalue)


# Prueba de chi-cuadrado

# Defino el umbral de longitud para dividir los datos en dos grupos
umbral_longitud = 59 # valor de la mediana

# Dividir los datos en dos grupos según la longitud
peptidos_pequenos = df[df['len'] <= umbral_longitud]
peptidos_grandes = df[df['len'] > umbral_longitud]

# Calcular las proporciones de "proteínas hipotéticas" en cada grupo
proporcion_pequenos = peptidos_pequenos[peptidos_pequenos['product'] == 'hypothetical protein'].shape[0] / peptidos_pequenos.shape[0]
proporcion_grandes = peptidos_grandes[peptidos_grandes['product'] == 'hypothetical protein'].shape[0] / peptidos_grandes.shape[0]

# Realizar una prueba de chi-cuadrado de independencia para comparar las proporciones
tabla_contingencia = pd.crosstab(df['len'] <= umbral_longitud, df['product'] == 'hypothetical protein')
_, p_valor, _, _ = chi2_contingency(tabla_contingencia)

# Imprimir los resultados
print("Proporción de 'proteínas hipotéticas' en péptidos más pequeños que la mediana:", proporcion_pequenos)
print("Proporción de 'proteínas hipotéticas' en péptidos más grandes que la mediana:", proporcion_grandes)
print("Valor p de la prueba de chi-cuadrado:", p_valor)

# Interpretación del resultado
nivel_significancia = 0.05 
if p_valor < nivel_significancia:
    print("Hay una diferencia significativa en las proporciones de 'proteínas hipotéticas' entre péptidos pequeños y grandes.")
else:
    print("No hay evidencia suficiente para afirmar que hay diferencias significativas en las proporciones de 'proteínas hipotéticas' entre péptidos pequeños y grandes.")


# ANOVA

conserved = df[df['Results_column'] == 'Conserved']['len'].dropna()
core = df[df['Results_column'] == 'Core']['len'].dropna()
soft_core = df[df['Results_column'] == 'Soft core']['len'].dropna()
shell = df[df['Results_column'] == 'Shell']['len'].dropna()
cloud = df[df['Results_column'] == 'Cloud']['len'].dropna()
probably_not_conserved = df[df['Results_column'] == 'Probably not conserved']['len'].dropna()


# Realizar el test ANOVA
anova_result = f_oneway(conserved, core, soft_core, shell, cloud, probably_not_conserved)

print("Resultado del ANOVA:")
print("Valor F:", anova_result.statistic)
print("Valor p:", anova_result.pvalue)

if anova_result.pvalue < 0.05:
    print("Hay diferencias significativas entre al menos dos grupos.")
else:
    print("No se encontraron diferencias significativas entre los grupos.")


sns.boxplot(x='Results_column', y='len', data=df)
plt.title('Distribución de Longitudes por Categoría de Conservación')
plt.xlabel('Categoría de Conservación')
plt.ylabel('Longitud')
plt.xticks(rotation=45)
plt.show()