# Trabajo final de la materia Biología Computacional
Estudiantes: Cecilia Rodriguez y Zapata Gonzalo

## Marco Teórico

Los genomas albergan una gran cantidad de Marcos de Lectura Abiertos (ORFs), sin embargo, determinar cuáles serán traducidos, especialmente en el caso de los ORFs pequeños que pueden ocurrir por casualidad, no es trivial, lo que complica el proceso de anotación de secuencias codificantes de proteínas. Las tecnologías modernas de secuenciación y proteómica han revelado que algunos ORFs pequeños (smORFs, por sus siglas en inglés) son activamente transcritos y traducidos en péptidos de menos de 70 aminoácidos (SEPs, SmORFs Encoded Peptids). Estos SEPs están involucrados en diversas funciones celulares, desde la esporulación y la división celular hasta el transporte, la regulación y la señalización, lo que los hace interesantes para su estudio. 

Aunque los avances recientes en los métodos de aprendizaje automático han logrado cierto progreso en la predicción de smORFs, aún queda mucho margen de mejora y el número real de smORFs en distintos organismos sigue siendo motivo de debate. Además, existe una falta significativa de caracterización de smORFs en varios grupos bacterianos, en particular en aquellos involucrados en interacciones con huéspedes eucariotas, como las bacterias promotoras del crecimiento de plantas en la rizosfera (PGPR), y bacterias patógenas tanto para plantas como para animales. En este estudio, realizamos un relevamiento de la información actual sobre smORFs ya anotados en los genomas secuenciados de *Sinorhizobium meliloti*, una alfa-proteobacteria del suelo capaz de fijar nitrógeno atmosférico cuando está en simbiosis con plantas leguminosas, que debido a su estilo de vida versátil y relevancia ecológica, sirve como un organismo modelo importante para estudiar la regulación genética en las alfaproteobacterias.

## Diseño y metodología del experimento

Para poder estudiar los smORFs de *S.meliloti*, se diseño un pipeline para descargar todos los genomas que se encuentren completos para esta bacteria en la base de datos RefSeq. Una vez descargados estos genomas, se creó un script de Python para poder extraer información de todas las secuencias codificantes de 70 residuos de largo o menos. Con esta información se creó un archivo csv que contiene para cada smORF presente el genoma correspondiente, el locus tag, la secuencia del péptido codificado, la longitud de esa secuencia, el producto anotado para ese péptido y la función en caso de estar disponible. Además, se agrego a la tabla los resultados de 3 programas distintos que se usaron para poder estudiar la conservación de estos genes a lo largo de los distintos genomas. Una vez agregados estos valores, se agregó una nueva columna con un resultado "consenso" a partir de los 3 programas. 

Con la tabla obtenida, se propone estudiar el tamaño promedio de estos genes, la información anotada y el grado de conservación de cada gen.

## 1) Descripción y Análisis de los datos

Se estudió el tamaño de los péptidos codificados para cada smORF. Esta variable solo puede tomar números enteros entre 0 y 70 por lo que se trata de una variable discreta. La forma elegida para representar el tamaño de estos péptidos fue a través de un histograma. 

![image](https://github.com/GonzaloMartinZapata/biologia_computacional/assets/173167223/af0bf620-901b-4681-8ea2-e5ef1df764e2)

En este gráfico se puede ver que la mayoría de los smORFs codifican para péptidos con secuencias más largas que 40 péptidos y que al considerar tamaños más grandes, mas genes son encontrados. De la misma forma, se puede decir que son muy pocos los genes que codifican péptidos de menos de 40 aminoácidos. No queda claro si esto se debe a que estos genes son poco frecuentes o a que los métodos de detección y anotación actuales no son capaces de reconocer estas pequeñas secuencias.

También, se propuso estudiar la distribución de estos genes por genoma. Para ello, se realizó un gráfico de tipo boxplot.

![image](https://github.com/GonzaloMartinZapata/biologia_computacional/assets/173167223/1d6c4511-785e-413a-ab9d-1ac10bae3047)

De este gráfico se observa una correlación con el anterior, al considerar tamaños más grandes, más smORFs son encontrados por genoma. A la vez, se puede decir que existe una varaicón en la cantidad de estos genes por genoma como se puede ver por el tamaño de las cajas.

Por último, se graficó la distribución del grado de conservación de estos genes en un gráfico de tipo torta, teniendo en cuenta los valores de la columna consenso.

![image](https://github.com/GonzaloMartinZapata/biologia_computacional/assets/173167223/6f02c527-4632-4a08-9128-89f0aa49df2a)

El siguiente paso fue calcular los parámetros de centralización y dispersión para el tamaño del producto de estos genes. Entre los parametros de centralización se calculo: la media, la mediana, la moda y el percentil 50. Los parametros de dispersión calculados fueron: el rango, el recorrido intercuartílico, la desviación típica y la varianza.

```python
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

Los resultados obtenidos fueron:

La media es 57.44264761188417
La mediana es 59.0
la moda es: ModeResult(mode=70, count=495)
El percentil 50.0 es 59.0
El rango es 56
El recorrido intercuartílico es 15.0
La desviación típica es 9.615593081095323
La varianza es 92.45963030120825

   ```

## 2) Asimetría y curtosis

Se calcularon también el coeficiente de asimetria de Fisher y el coeficiente de curtosis.

 coeficiente de asimetría:
  ```python
skewness = df["len"].skew(axis=0, skipna=True, numeric_only=False)
 ```
El valor obtenido fue de -0.9114654884958558. Al ser neagtivo se puede decir que la distribución de tamaños de smORFs en S.meliloti tiene una asimetría negativa (hacia la izquierda)

coeficiente de curtosis:
  ```python
kurtosis = df["len"].kurt(axis=0, skipna=True, numeric_only=False)
 ```
En este caso el resultado fue 1.0995996120417915 lo que indica una distribución leptocúrtica. Esto indica una mayor concentración de valores cerca de la media.

## 3) Intervalos de confianza

Se determinaron los intervalos de confianza para el promedio de smORFs presentes por genoma en total y también para aquellos que puedan ser considerados como "probablemente conservados". Se consideran probablemente conservados a aquellos que fueron clasificados en las categorías "conserved", "core" y "soft core" en la columna de resultados. Estos genes estarían presentes en el 100% de los genomas estudiados, entre 99 y 95% o entre 50 y 95% respectivamente.

Para este cálculo se uso la variable t de student ya que el n era 28 (el total de genomas estudiados). El nivel de confianza utilizado fue 0.95 y los grados de libertad fueron 27 (n-1). Los resultados obtenidos fueron:

-La media del promedio de smORFS por genoma es: 284.89285714285717 <br>
-La mediana del promedio de smORFs por genoma es: 277.0 <br>
-La desviación típica del promedio de smORFS por genoma es: 32.80169291669346 <br>
-La varianza del promedio de smORFs por genoma es: 1075.9510582010582 <br>
-Intervalo de confianza al 95.0%: (272.173688444495, 297.61202584121935) <br>
-La media de smORFs conservados por genoma es: 188.5 <br>
-La mediana del promedio de smORFs conservados por genoma es: 192.0 <br>
-La desviación típica del promedio de smORFS conservados por genoma es: 11.516493807996843 <br>
-La varianza del promedio de smORFs conservados por genoma es: 132.62962962962962 <br>
-Intervalo de confianza al 95.0% para smORFs conservados: (184.0343698897014, 192.9656301102986) <br>

## 4) Test de Hipótesis

Se propuso estudiar como están anotadas estos SEPs. Para esto se realizó un gráfico de tipo torta mostrando los distintos productos codificados por estos genes.

![image](https://github.com/GonzaloMartinZapata/biologia_computacional/assets/173167223/25e0f2c6-0df1-4448-b126-307f986b07a0)

Se puede observar que aproximadamente dos tercios de estos péptidos están anotados como proteínas hipotéticas, lo que significa que no existe validación experimental que compruebe que estas proteínas se expresan. Entre aquellos SEPs que si cuentan con una anotación se encuentra una gran variedad de dominios y familias distintas que muestran los diversos roles que pueden llevar acabo estas pequeñas proteínas.

A continuación se planteao un test de hipótesis para comprobar si existe una diferencia en la anotación de estas proteínas dada por el tamaño. Las técnicas de laboratorio usadas generalmente tienden a no ser específicas para proteínas tan chicas, esto sumado a la falta de caracterización que dificulta las estimaciones por métodos bioinformáticos, nos llevó a plantear si esta dificultad se hace más notoria con péptidos más cortos. Para este análisis se plantearon las siguientes hipótesis:

H0: No hay diferencia en las proporciones de "proteínas hipotéticas" entre los péptidos pequeños (longitud menor o igual a la mediana) y los péptidos grandes (longitud mayor a la mediana).

H1: Hay una diferencia en las proporciones de "proteínas hipotéticas" entre los péptidos pequeños y los péptidos grandes.

Se usó en este caso la prueba de chi cuadrado. Los resultados obtenidos fueron: 

Proporción de 'proteínas hipotéticas' en péptidos más pequeños que la mediana: 0.7482551143200963 <br>
Proporción de 'proteínas hipotéticas' en péptidos más grandes que la mediana: 0.5722135007849294 <br>
Valor p de la prueba de chi-cuadrado: 6.474964936751062e-62 <br>
Hay una diferencia significativa en las proporciones de 'proteínas hipotéticas' entre péptidos pequeños y grandes. <br>

Lo que confirma que para péptidos más pequeños la falta de información y caracterización es mayor.

Por otro lado, se realizó un test de ANOVA para comparar si existen diferencias significativas en los tamaños de los SEPs teniendo en cuenta su grado de conservación. En este caso, los resultados fueron:

Valor F: 22.6498951969155 <br>
Valor p: 1.2290935614758286e-22 <br>
Hay diferencias significativas entre al menos dos grupos <br>

Para representar estos resultados de forma más clara, se realizó un nuevo boxplot con estos resultados:

![image](https://github.com/GonzaloMartinZapata/biologia_computacional/assets/173167223/791881cf-02c2-404e-a936-6cdd9c6d627b)

