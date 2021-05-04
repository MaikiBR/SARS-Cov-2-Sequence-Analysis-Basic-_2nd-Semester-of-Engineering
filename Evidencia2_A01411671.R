# Evidencia 2 | Proyecto integrador
# Análisis de biología computacional (Gpo 519)
# Miguel Ángel Bermea Rodríguez | A01411671

## PARTE 1 - Video (Esta en el comentario de Canvas)
# https://drive.google.com/drive/folders/1shCcmLZwI3JlVDujccZKc3gwWmfnNGnU?usp=sharing
# https://youtu.be/i9Tod5911eg

## PARTE 2 - Código
cat('Evidencia 2 | Proyecto integrador\n','Análisis de biología computacional (Gpo 519)\n',
    'Miguel Ángel Bermea Rodríguez | A01411671\n','PARTE 2 - Código')

# 1. Secuencias
cat('1. Secuencias')
# US - MZ058005.1 |Severe acute respiratory syndrome coronavirus 2 isolate 
# SARS-CoV-2/human/USA/TX-CDC-STM-000053264/2021, complete genome

# India - MW881790.1 |Severe acute respiratory syndrome coronavirus 2 isolate 
# SARS-CoV-2/human/IND/THSTI_UK_18022021_CSIR-IGIB/2021, complete genome

# Brasil - MW592707.1 |Severe acute respiratory syndrome coronavirus 2 isolate 
# SARS-CoV-2/human/BRA/MASP2C844R2/2020, complete genome

# Francia - MW580244.1 |Severe acute respiratory syndrome coronavirus 2 isolate 
# SARS-CoV-2/human/FRA/ZA-1/2021, complete genome

# Turquia - MW308549.1 |Severe acute respiratory syndrome coronavirus 2 isolate 
# SARS-CoV-2/human/TUR/DENIZLI_PAU/2020, complete genome

cat('# US - MZ058005.1 |Severe acute respiratory syndrome coronavirus 2 isolate 
# SARS-CoV-2/human/USA/TX-CDC-STM-000053264/2021, complete genome

# India - MW881790.1 |Severe acute respiratory syndrome coronavirus 2 isolate 
# SARS-CoV-2/human/IND/THSTI_UK_18022021_CSIR-IGIB/2021, complete genome

# Brasil - MW592707.1 |Severe acute respiratory syndrome coronavirus 2 isolate 
# SARS-CoV-2/human/BRA/MASP2C844R2/2020, complete genome

# Francia - MW580244.1 |Severe acute respiratory syndrome coronavirus 2 isolate 
# SARS-CoV-2/human/FRA/ZA-1/2021, complete genome

# Turquia - MW308549.1 |Severe acute respiratory syndrome coronavirus 2 isolate 
# SARS-CoV-2/human/TUR/DENIZLI_PAU/2020, complete genome')

# 2. Calcula la longitud de las secuencias que incluyas.
cat('2. Calcula la longitud de las secuencias que incluyas.')

library(seqinr)

cat('Secuencia US | MZ058005')
secuenciaUS<-read.fasta(file.choose()) #Archivo: "US_TX.fasta"
US_seq<-secuenciaUS[[1]]
cat('Longitud de la Secuencia US | MZ058005:\n',length(US_seq)) # Longitud de toda la secuencia
cat('Magnitud de cada base de la Secuencia US | MZ058005:\n')
table(US_seq) # Magnitud de cada base

cat('Secuencia India | MW881790')
secuenciaIndia<-read.fasta(file.choose()) #Archivo: "India.fasta"
India_seq<-secuenciaIndia[[1]]
cat('Longitud de la Secuencia India | MW881790:\n',length(India_seq)) # Longitud de toda la secuencia
cat('Magnitud de cada base de la Secuencia India | MW881790:\n')
table(India_seq) # Magnitud de cada base

cat('Sequencia Brasil | MW592707')
secuenciaBrasil<-read.fasta(file.choose()) #Archivo: "Brazil.fasta"
Brasil_seq<-secuenciaBrasil[[1]]
cat('Longitud de la Secuencia Brasil | MW592707:\n',length(Brasil_seq)) # Longitud de toda la secuencia
cat('Magnitud de cada base de la Secuencia Brasil | MW592707:\n')
table(Brasil_seq) # Magnitud de cada base

cat('Sequencia Francia | MW580244')
secuenciaFrancia<-read.fasta(file.choose()) #Archivo: "Francia.fasta"
Francia_seq<-secuenciaFrancia[[1]]
cat('Longitud de la Secuencia Francia | MW580244:\n',length(Francia_seq)) # Longitud de toda la secuencia
cat('Magnitud de cada base de la Secuencia Francia | MW580244:\n')
table(Francia_seq) # Magnitud de cada base

cat('Sequencia Turquia | MW308549')
secuenciaTurquia<-read.fasta(file.choose()) #Archivo: "Turquia.fasta"
Turquia_seq<-secuenciaTurquia[[1]]
cat('Longitud de la Secuencia Turquia | MW308549:\n',length(Turquia_seq)) # Longitud de toda la secuencia
cat('Magnitud de cada base de la Secuencia Turquia | MW308549:\n')
table(Turquia_seq) # Magnitud de cada base

# 3. Crea una sola gráfica donde se comparen el número de bases de ADN que 
# componen todas las variantes del virus.
cat('3. Crea una sola gráfica donde se comparen el número de bases de ADN que 
    componen todas las variantes del virus.')

 # Grafica comparacion de todas las secuencias

library(ggplot2)
Country = c(rep("Secuencia US",4),rep("Secuencia India",4),rep("Sequencia Brasil",4),rep("Sequencia Francia",4),rep("Sequencia Turquia",4))
Nucleotides = rep(c("Adenine","Cytosine","Guanine","Thymine"),5)
Bases = c(8908,5475,5850,9583,8912,5477,5849,9593,8923,5485,5856,9598,8900,5473,5850,9589,8899,5474,5852,9588)
data = data.frame(Country,Nucleotides,Bases)
ggplot(data, aes(fill=Nucleotides, y=Bases, x=Country)) + 
  geom_bar(position="stack", stat="identity")
  
  # Graficas individuales

grafica<-function(seq1,seq2,seq3,seq4,seq5){
  par(mfrow=c(2,3))
  barplot(table(seq1),col=1:4,main="US")
  barplot(table(seq2),col=1:4,main="India")
  barplot(table(seq3),col=1:4,main="Brasil")
  barplot(table(seq4),col=1:4,main="Francia")
  barplot(table(seq5),col=1:4,main="Turquia")
}

grafica(US_seq,India_seq,Brasil_seq,Francia_seq,Turquia_seq)


# 4. Agrega un análisis jerárquico global (árbol filogenetico [las 4 formas]) 
# obtenido de las secuencias que se seleccionaron para estudiar.
cat('4. Agrega un análisis jerárquico global (árbol filogenetico [las 4 formas]) 
    obtenido de las secuencias que se seleccionaron para estudiar.')

  # Tipo 1
dev.off()
library(ape)
text.string<-
  "((MW881790|India,MZ058005|US),MW592707|Brasil,(MW308549|Turquia,MW580244|Francia));"
vert.tree<-read.tree(text=text.string)
plot(vert.tree,no.margin=TRUE,edge.width = 2)
  # Tipo 2
library(phytools)
roundPhylogram(vert.tree)
  # Tipo 3
plot(unroot(vert.tree),type="unrooted",no.margin=TRUE,lab4ut="axial",
     edge.width=2)
  # Tipo 4
vert.tree
str(vert.tree)

tree<-read.tree(text="((MW881790|India,MZ058005|US),MW592707|Brasil,(MW308549|Turquia,MW580244|Francia));")
plotTree(tree,offset=1)
tiplabels()
nodelabels()


# 5. Interpretación escrita de tus gráficas y tus conclusiones según el caso de 
# estudio que seleccionaste. No olvides sustentar tus argumentos con las 
# lecturas que realizaste.

# Analizar las secuencias de SARS-CoV-2 reportadas en los 20 países con más casos 
# reportados. Y puedes tratar de responder a la pregunta: 
# ¿Son muy diferentes las variantes entre cada país? ¿Es diferente el SARS-CoV-2 
# entre las diferentes poblaciones: Asiática, Hispana, Europea o Africana?
cat('5. Interpretación escrita de tus gráficas y tus conclusiones según el caso de 
estudio que seleccionaste. No olvides sustentar tus argumentos con las 
    lecturas que realizaste.')

cat('En este trabajo se decidio investigar el primer caso en el que investigue 
secuencia de SARS-CoV-2 reportadas en los 5 paises con mas casos reportados 
(Secuencias obtenidas de NCBI). Primeramente y respondiendo a la primera 
interrogante, las diferencias entre a variente de cada pais no es mucha, 
esto se puede apreciar primero por la longitud de las secuencia, ya que las 5 
secuencias estan en entorno a 29800 (especificamente un rango de [29812 a 29862]). 
Ademas en todas las secuencias la Adenina y la Timina son las de mayor frecuencia.')

cat('Apoyando esta respuesta tambien esta el histograma (sobretodo el tipo 4 que nos
muestra los nodos de una manera visual), ya que el que secuencias compartan
nodo quiere decir que se trata de una pequeña mutacion o variante, por ejemplo,
en este caso la secuencia de Francia y la de Turquia comparten nodo, al igual que
la secuencia de US e India, por lo tanto estas son bastante similares entre si.')

cat('Y la grafica comparativa de las secuencias confirma esto, ya que podemos apreciar
lo que mencione anteriormente, la mayor frecuencia es de Adenina y Timina, y
ademas podemos apreciar visualmente la gran similitud que existen entre esas variantes
del Sars-CoV-2 en los 5 paises con mas casos reportados.')

cat('Por lo que extrapolando lo contestado en la primera pregunta y respondiedo ahora la
segunda interrogante, el SARS-CoV-2 si es ligeramente diferente entre las
poblaciones, mas se tratan de diferencias minimas (no son muy diferentes), 
y son mas las similitudes que diferencias entre las diferentes poblaciones, 
por lo menos, de estos 5 paises que son los que mas casos han reportado, 
lo que sugiere una tendencia de esta variante de ser la que mas casos reporte,
ademas todas conclusiones se ven completamente sustentadas con lo mencionado 
en los articulos leidos (PUBMED: Evolucion del coronavirus y su compartamiento
en y desarrollo en distintos paises y entornos) tanto para la PARTE 1 como para esta PARTE 2 de la Evidencia
del Proyecto Integrador.')


# 6. Referencias de la PARTE 1 y de la PARTE 2
cat('6. Referencias de la PARTE 1 y de la PARTE 2')

'Wu, F., Zhao, S., Yu, B., Chen, Y. M., Wang, W., Song, Z. G., Hu, Y., Tao, Z. 
W., Tian, J. H., Pei, Y. Y., Yuan, M. L., Zhang, Y. L., Dai, F. H., Liu, Y., 
Wang, Q. M., Zheng, J. J., Xu, L., Holmes, E. C., & Zhang, Y. Z. (2020). 
A new coronavirus associated with human respiratory disease in China. Nature, 
579(7798), 265-269. https://doi.org/10.1038/s41586-020-2008-3'

'Rabaan, A. A., Al-Ahmed, S. H., Haque, S., Sah, R., Tiwari, R., Malik, Y. S., 
Dhama, K., Yatoo, M. I., Bonilla-Aldana, D. K., & Rodriguez-Morales, A. J. (2020). 
SARS-CoV-2, SARS-CoV, and MERS-COV: A comparative overview. Le infezioni in medicina, 28(2), 174-184.'

'Li, C., Yang, Y., & Ren, L. (2020). Genetic evolution analysis of 2019 novel 
coronavirus and coronavirus from other species. Infection, genetics and evolution : 
journal of molecular epidemiology and evolutionary genetics in infectious diseases, 
82, 104285. https://doi.org/10.1016/j.meegid.2020.104285'