Este repositorio contiene los desarrollos en el campo de radio ocultación GNSS para la tesina de grado del alumno SAGGESE, ARIAN. 


Descripción por carpeta

PROPAGADOR SIMPLE

Contiene los códigos para simular cualquier órbita dado un archivo YUMA o por descripción manual, también se contienen los script de captación y clasificación de eventos de RO

LAZOS DE SEGUIMIENTO 

Contiene los códigos de lazos de seguimiento sin filtro de Kalman, utilizan PLL

KALMAN

Tiene algunos scripts de implementación del filtro de Kalman para escenarios particulares, se utilizaban para comprender el funcionamiento

GENERADOR DE SEÑAL RO

Acá se encuentran los scripts de generación de señal de RO (para luego generar una señal GPS), este script *GenModRO* permite generar señales GPS con características configurables como la frecuencia de muestreo si se quiere o no ensanchar por código, etc. Para luego poder utilizar los códigos de seguimiento por filtro de Kalman. 
