### Precision

homology <- matrix(c(1,0.263,5,0.307,10,0.330,20,0.350,30,0.384,40,0.398,50,0.400),ncol=2,byrow=TRUE)
colnames(homology) <- c("NAMP/AMP","Precision")
rownames(homology) <- NULL
homology <- as.table(homology)

amp_scan <- matrix(c(1,0.923,5,0.953,10,0.980,20,0.992,30,0.994,40,0.986,50,0.989),ncol=2,byrow=TRUE)
colnames(amp_scan) <- c("NAMP/AMP","Precision")
rownames(amp_scan) <- NULL
amp_scan <- as.table(amp_scan)

macrel <- matrix(c(1,0.949,5,0.983,10,0.994,20,0.994,30,1.000,40,1.000,50,1.000),ncol=2,byrow=TRUE)
colnames(macrel) <- c("NAMP/AMP","Precision")
rownames(macrel) <- NULL
macrel <- as.table(macrel)

iamp2L <- matrix(c(1,0.917,5,0.963,10,0.974,20,0.982,30,0.989,40,0.991,50,0.994),ncol=2,byrow=TRUE)
colnames(iamp2L) <- c("NAMP/AMP","Precision")
rownames(iamp2L) <- NULL
iamp2L <- as.table(iamp2L)

svg(filename="Pr_plot.svg", 
    width=5, 
    height=4, 
    pointsize=12)
plot(macrel,type = "o",col = "blue", xlab = "non-AMP:AMP", ylab = "Precision", ylim=c(0.9,1.01))
lines(iamp2L, type="o", col="green")
lines(amp_scan, type="o", col="orange")
dev.off()

## Accuracy

homology <- matrix(c(1,0.222,5,0.318,10,0.352,20,0.384,30,0.421,40,0.432,50,0.437),ncol=2,byrow=TRUE)
colnames(homology) <- c("NAMP/AMP","Precision")
rownames(homology) <- NULL
homology <- as.table(homology)

amp_scan <- matrix(c(1,0.908,5,0.923,10,0.894,20,0.889,30,0.845,40,0.856,50,0.850),ncol=2,byrow=TRUE)
colnames(amp_scan) <- c("NAMP/AMP","Precision")
rownames(amp_scan) <- NULL
amp_scan <- as.table(amp_scan)

macrel <- matrix(c(1,0.926,5,0.893,10,0.849,20,0.804,30,0.763,40,0.734,50,0.718),ncol=2,byrow=TRUE)
colnames(macrel) <- c("NAMP/AMP","Precision")
rownames(macrel) <- NULL
macrel <- as.table(macrel)

iamp2L <- matrix(c(1,0.931,5,0.922,10,0.907,20,0.870,30,0.860,40,0.841,50,0.829),ncol=2,byrow=TRUE)
colnames(iamp2L) <- c("NAMP/AMP","Precision")
rownames(iamp2L) <- NULL
iamp2L <- as.table(iamp2L)

svg(filename="Acc_plot.svg", 
    width=5, 
    height=4, 
    pointsize=12)
plot(homo,type = "o",col = "red", xlab = "non-AMP:AMP", ylab = "Accuracy", ylim=c(0,1))
lines(macrel, type="o", col="blue")
lines(iamp2L, type="o", col="green")
lines(amp_scan, type="o", col="orange")
legend(40, 0.25, legend = c("Homology", "Macrel", "AMP Scanner", "iAMP-2L"), col=c("red","blue","green","orange"), pch=1, cex=0.5)
dev.off()





