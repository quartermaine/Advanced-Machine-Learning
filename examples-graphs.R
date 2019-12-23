library(bnlearn)
i
dag<-empty.graph(nodes=c("A","S","E","O","R","T"))
dag

dag<-set.arc(dag,from="A",to="E")
dag<-set.arc(dag,from = "S",to="E")

dag<-set.arc(dag,from="E",to="O")
dag<-set.arc(dag,from="E",to="R")


dag<-set.arc(dag,from="O",to="T")
dag<-set.arc(dag,from="R",to="T")

dag

modelstring(dag)

nodes(dag)
arcs(dag)


# Dag matrix insertion ----------------------------------------------------

dag2 <- empty.graph(nodes = c("A", "S", "E", "O", "R", "T"))
arc.set <- matrix(c("A", "E",
                    "S", "E",
                    "E", "O",
                    "E", "R",
                    "O", "T",
                    "R", "T"),byrow = TRUE, ncol = 2,
                  dimnames = list(NULL, c("from", "to")))
arcs(dag2) <- arc.set

dag2

all.equal(dag,dag2)

try(set.arc(dag,from="T",to="E"))


A.lv<-c("young", "adult", "old")
S.lv<-c("M", "F")
E.lv<-c("high", "uni")
O.lv<-c("emp", "self")
R.lv<-c("small", "big")
T.lv<-c("car", "train", "other")

A.prob <- array(c(0.30, 0.50, 0.20), dim = 3,
                dimnames = list(A = A.lv))
A.prob

S.prob <- array(c(0.60, 0.40), dim = 2,
                dimnames = list(S = S.lv))
S.prob

O.prob <- array(c(0.96, 0.04,0.92,0.08),dim=c(2,2),
                dimnames = list(O=O.lv,E=E.lv))

O.prob                  
                  
R.prob<-array(c(0.25,0.75,0.20,0.80),dim=c(2,2),
              dimnames=list(R=R.lv,E=E.lv))
  
R.prob  

E.prob <- array(c(0.75, 0.25, 0.72, 0.28, 0.88, 0.12, 0.64,
                  0.36, 0.70, 0.30, 0.90, 0.10), dim = c(2, 3, 2),
                dimnames = list(E = E.lv, A = A.lv, S = S.lv))
E.prob

T.prob <- array(c(0.48, 0.42, 0.10, 0.56, 0.36, 0.08, 0.58,
                  0.24, 0.18, 0.70, 0.21, 0.09), dim = c(3, 2, 2),
                dimnames = list(T = T.lv, O = O.lv, R = R.lv))
T.prob

dag3 <- model2network("[A][S][E|A:S][O|E][R|E][T|O:R]")
all.equal(dag, dag3)

cpt <- list(A = A.prob, S = S.prob, E = E.prob, O = O.prob,
              R = R.prob, T = T.prob)
bn <- custom.fit(dag, cpt)

nparams(bn)
arcs(bn)


bn$R

R.cpt<-coef(bn$R)
bn


survey <- read.table("survey.txt", header = TRUE)
head(survey)
bn.mle<-bn.fit(dag,data=survey,method="mle")

prop.table(table(survey[, c("O", "E")]), margin = 2)

bn.mle$O


bn.bayes <- bn.fit(dag, data = survey, method = "bayes",iss = 10)

bn.bayes$O

ci.test("T", "E", c("O", "R"), test = "mi", data = survey)

ci.test("T", "E", c("O", "R"), test = "x2", data = survey) # Pearson's X^2

arc.strength(dag, data = survey, criterion = "x2")

score(dag, data = survey, type = "bic") # BIC
score(dag, data = survey, type = "bde", iss = 10) # BDE

dag4 <- set.arc(dag, from = "E", to = "T")
nparams(dag4, survey)
score(dag4, data = survey, type = "bic")

rnd <- random.graph(nodes = c("A", "S", "E", "O", "R", "T"))
modelstring(rnd)
score(rnd, data = survey, type = "bic")

learned<-hc(survey)
modelstring(learned)
score(learned,data=survey,type="bic")

arc.strength(learned, data = survey, criterion = "bic")
arc.strength(dag, data = survey, criterion = "bic")

dsep(dag, x = "S", y = "R")
dsep(dag, x = "O", y = "R")

path(dag,from="S",to="R")
plot(dag)

dsep(dag, x = "S", y = "R", z = "E")


library(gRain)

junction <- compile(as.grain(bn))

querygrain(junction,nodes="T")$T
jsex<-setEvidence(junction,nodes="S",states="F")
querygrain(jsex,nodes="T")$T

jres<-setEvidence(junction,nodes="R",states="small")
querygrain(jres,nodes="T")$T


jedu <- setEvidence(junction, nodes = "E", states = "high")
SxT.cpt <- querygrain(jedu, nodes = c("S", "T"),
                      type = "joint")
SxT.cpt

querygrain(jedu, nodes = c("S", "T"), type = "conditional")
dsep(bn, x = "S", y = "T", z = "E")



SxT.ct = SxT.cpt * nrow(survey)
chisq.test(SxT.ct)


cpquery(bn, event = (S == "M") & (T == "car"),
        evidence = (E == "high"))

cpquery(bn, event = (S == "M") & (T == "car"),
        evidence = (E == "high"),n=10^6)


cpquery(bn, event = (S == "M") & (T == "car"),
        evidence = list(E = "high"), method = "lw")


SxT <- cpdist(bn, nodes = c("S", "T"),
              evidence = (E == "high"))
head(SxT)

graphviz.plot(dag)


hlight <- list(nodes = nodes(dag), arcs = arcs(dag),
               col = "grey", textCol = "grey")

pp <- graphviz.plot(dag, highlight = hlight)
edgeRenderInfo(pp) <-list(col = c("S~E" = "black", "E~R" = "black"),
                          lwd = c("S~E" = 3, "E~R" = 3))

nodeRenderInfo(pp) <-list(col = c("S" = "black", "E" = "black", "R" = "black"),
                          textCol = c("S" = "black", "E" = "black", "R" = "black"),
                          fill = c("E" = "grey"))
renderGraph(pp)

bn.fit.barchart(bn.mle$T, main = "Travel",
                xlab = "Pr(T | R,O)", ylab = "")


Evidence<-factor(c(rep("Unconditional",3), rep("Female", 3),
                    rep("Small City",3)),
                  levels = c("Unconditional", "Female", "Small City"))
Travel<-factor(rep(c("car", "train", "other"), 3),
               levels = c("other", "train", "car"))
distr<-data.frame(Evidence = Evidence, Travel = Travel,
                  Prob = c(0.5618, 0.2808, 0.15730, 0.5620, 0.2806,
                           0.1573, 0.4838, 0.4170, 0.0990))

barchart(Travel ~ Prob | Evidence, data = distr,
         layout = c(3, 1), xlab = "probability",
         scales = list(alternating = 1, tck = c(1, 0)),
         strip = strip.custom(factor.levels =c(expression(Pr(T)),expression(Pr({T} * " | " * {S == F})),
                                  expression(Pr({T} * " | " * {R == small})))),
         panel = function(...) {
           panel.barchart(...)
           panel.grid(h = 0, v = -1)})
  
