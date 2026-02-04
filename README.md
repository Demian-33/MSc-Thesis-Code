En este repositorio se encuentran dos contenidos principales:



1. En la carpeta `Codigos`, se encuentran scripts para reproducir los análisis, esto es, ajustes y estimaciones basadas en la distribución \*a posteriori\*.



Dentro de `Codigos`, también se encuentran scripts para reproducir diagramas y las gráficas básicas, por ejemplo: densidad normal asimétrica, divergencia KL estimada, diagramas.



2\. En la carpeta `LaTeX` se dispone de código `.tex` para construir un documento en formato `.pdf`. El documento raíz es `Main-Doc.tex`. Es necesario modificar algunas rutas para evitar problemas en la compilación:



Por ejemplo, reemplazar



`\\includegraphics{../Figuras/c-ii/SN\_density\_pos.pdf}`



en lugar de



`\\includegraphics{Figuras/c-ii/SN\_density\_pos.pdf}`

