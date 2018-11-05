2
=

rna\_counts &lt;- A3

dim(rna\_counts)

mean\_fun &lt;- function (x=rna\_counts, y=FALSE) { if (y == TRUE) { x
&lt;- (log2(x+0.00000001))} return(mean(x))}

mean\_fun(rna\_counts\[,2\])  
mean\_fun(rna\_counts\[,2\], y = TRUE)

3
=

dim\_data &lt;- ncol(rna\_counts) stored\_means &lt;- rep(0,
dim\_data-1)

for (i in 1:dim\_data) { stored\_means\[i-1\] &lt;-
mean\_fun(rna\_counts\[,i\], y=FALSE) print(colnames(rna\_counts\[i\]))
print(stored\_means\[i-1\]) }

4
=

sapply(rna\_counts\[,2:56\], mean\_fun)

TimeTest2 &lt;- sapply(rna\_counts\[,2:56\], mean\_fun)

TimeTest1 &lt;- for (i in 1:dim\_data) { stored\_means\[i-1\] &lt;-
mean\_fun(rna\_counts\[,i\], y=FALSE) print(colnames(rna\_counts\[i\]))
print(stored\_means\[i-1\]) }

library(microbenchmark)

microbenchmark(TimeTest2) microbenchmark(TimeTest1)

The sapply function (TimeTest2 = Question 4) ran faster than the for loop (TimeTest1 = Question 3)
==================================================================================================

5
=

colMeans(rna\_counts\[,2:56\])

6
=

nrow(rna\_counts) stored\_gene\_means &lt;-
rowMeans(rna\_counts\[,2:56\]) print(stored\_gene\_means)

7
=

We are very interested in what is going on in the head horns between small males and large males.
=================================================================================================

Using the type of tools you have written (feel free to modify as you need, but show the new functions) calculate the mean expression for the subset of columns for large and small male head horns.
===================================================================================================================================================================================================

Note you are calculating means on a gene by gene basis, NOT sample by sample.
=============================================================================

mean\_gene\_diff &lt;- function (x) { smM &lt;- subset(x, select =
(grepl("sm\_male\_hdhorn*", names(x)))) smM\_mean &lt;- rowMeans(smM)
lgM &lt;- subset(x, select = (grepl("lg\_male\_hdhorn*", names(x))))
lgM\_mean &lt;- rowMeans(lgM)

Now calculate the mean difference (again gene by gene) between large male and small males (for head horns).
===========================================================================================================

i.e. first calculate the mean expression among individuals who are large males (head horns), ditto for the small males, and calculate their difference.
=======================================================================================================================================================

print(lgM\_mean - smM\_mean) }

rna\_gene\_diff &lt;- mean\_gene\_diff(rna\_counts)

8
=

Using the basic plot function (although you can use ggplot2 if you prefer), plot the mean expression of each gene on the X axis, and the difference in expression values on the Y axis.
=======================================================================================================================================================================================

library(ggplot2) theme\_update(plot.title = element\_text(hjust = 0.5))

ggplot(rna\_counts, aes(x=stored\_gene\_means, y=rna\_gene\_diff)) +
geom\_point() + ggtitle("Mean Gene Expression vs. Mean Gene Expression
Diffrence of lgM and smM") + xlab("Mean RNA Gene Expression") +
ylab("Mean RNA Expression Difference")

Now repeat, but with log2 transformed data. This is the basic idea of a MAplot.
===============================================================================

ggplot(rna\_counts, aes(x=stored\_gene\_means, y=rna\_gene\_diff)) +
geom\_point() + scale\_x\_continuous(trans='log2') +
scale\_y\_continuous(trans='log2') + ggtitle("Log of Mean Gene
Expression vs. Mean Gene Expression Diffrence of lgM and smM") +
xlab("Log Mean RNA Gene Expression") + ylab("Mean RNA Expression
Difference")
