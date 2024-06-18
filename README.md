# Trabajo final de la materia Biología Computacional
Estudiantes: Zapata Gonzalo

Marco Teórico

Los genomas albergan una gran cantidad de Marcos de Lectura Abiertos (ORFs), sin embargo, determinar cuáles serán traducidos, especialmente en el caso de los ORFs pequeños que pueden ocurrir por casualidad, no es trivial, lo que complica el proceso de anotación de secuencias codificantes de proteínas. Las tecnologías modernas de secuenciación y proteómica han revelado que algunos ORFs pequeños (smORFs, por sus siglas en inglés) son activamente transcritos y traducidos en péptidos de menos de 70 aminoácidos (SEPs, SmORFs Encoded Peptids). Estos SEPs están involucrados en diversas funciones celulares, desde la esporulación y la división celular hasta el transporte, la regulación y la señalización, lo que los hace interesantes para su estudio. 

Aunque los avances recientes en los métodos de aprendizaje automático han logrado cierto progreso en la predicción de smORFs, aún queda mucho margen de mejora y el número real de smORFs en distintos organismos sigue siendo motivo de debate. Además, existe una falta significativa de caracterización de smORFs en varios grupos bacterianos, en particular en aquellos involucrados en interacciones con huéspedes eucariotas, como las bacterias promotoras del crecimiento de plantas en la rizosfera (PGPR), y bacterias patógenas tanto para plantas como para animales. En este estudio, realizamos un relevamiento de la información actual sobre smORFs ya anotados en los genomas secuenciados de Sinorhizobium meliloti, una alfa-proteobacteria del suelo capaz de fijar nitrógeno atmosférico cuando está en simbiosis con plantas leguminosas, que debido a su estilo de vida versátil y relevancia ecológica, sirve como un organismo modelo importante para estudiar la regulación genética en las alfaproteobacterias.

Diseño y metodología del experimento

Para poder estudiar los smORFs de S.meliloti, se diseño un pipeline para descargar todos los genomas que se encuentren completos para esta bacteria en la base de datos RefSeq. Una vez descargados estos genomas, se creó un script de Python para poder extraer información de todas las secuencias codificantes que codifican para péptidos de 70 aminoácidos de largo o menos. Con esta información se creó un archivo csv que contiene para cada smORF presente el genoma correspondiente, el locus tag, la secuencia del péptido codificado, la longitud de esa secuencia, el producto anotado para ese péptido y la función en caso de estar disponible. Además, se agrego a la tabla, utilizando un script distinto, los resultados de 3 programas distintos que se usaron para poder estudiar la conservación de estos genes a lo largo de los distintos genomas. Una vez agregados estos valores, el mismo script agrega una nueva columna con un resultado "consenso" a partir de los 3 programas. 

Con la tabla obtenida, se propone estudiar el tamaño promedio de estos genes, la información anotada y el grado de conservación de cada gen.

1) Descripción y Análisis de los datos

Se estudió el tamaño de los péptidos codificados para cada smORF. Al estar analizando el tamaño del péptido codificado, esta variable solo puede tomar números enteros entre 0 y 70 por lo que se trata de una variable discreta. La forma elegida para representar el tamaño de estos péptidos fue a través de un histograma. 

