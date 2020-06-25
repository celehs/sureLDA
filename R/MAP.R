#' Internal function to do inverse logit transformation
#' @noRd
g.logit = function(xx){exp(xx)/(exp(xx)+1)}

#' Internal function to do logit transformation
#' @noRd
logit = function(xx){log(xx/(1-xx))}

#' Internal function to get the probabilities that match the prevalence
#' @noRd
Prob.Scale = function(pp, prev){
    logit_pp = logit(pp);
    logit_pp[logit_pp == -Inf] = min(logit_pp[logit_pp > -Inf]);
    logit_pp[logit_pp ==  Inf] = max(logit_pp[logit_pp < Inf])
    cc = uniroot(function(xx){mean(g.logit(logit_pp-xx))-prev},interval=c(-10,10))$root
    g.logit(logit_pp-cc)
}

#' Internal function to do MAP algorithm
#' @noRd
MAP_PheWAS_JS = function(dat=NULL, vname=NULL, yes.con=FALSE){

    vname.log = paste(vname,"_log",sep="")
    avgcount = apply(dat[,vname.log,drop=FALSE],1,mean)
    prior.all = NULL
    post.all = NULL
    group.all = NULL

    family = c( rep("poisson", length(vname)),rep("gaussian", length(vname.log)) )

    name.all = c(vname,vname.log)

    for(i in seq_along(name.all)){
        tmpfm = as.formula(paste(name.all[i],"~1"))
        if(yes.con){
            ii = i%%(length(vname))
            if(ii==0){ii = length(vname)}
            vname.log.out = setdiff(vname.log, vname.log[ii])
            tmpfm2 = as.formula(paste("~", paste(c("note",vname.log.out),collapse="+")))
            tmpfm2 = FLXPmultinom(tmpfm2)
        }else{
            tmpfm2 = FLXPconstant()
        }
        set.seed(1)
        n.clust = 1
        iter = 0
        while(n.clust < 2 & iter < 5){
            tmpfit = flexmix(tmpfm, k = 2,
                             model = FLXMRglmfix(fixed =~note,varFix=FALSE, family=family[i]),
                             concomitant=tmpfm2,control=list(tol = 1e-8), data=dat)
            n.clust = length(unique(tmpfit@cluster))
            iter = iter+1
        }

        ##if(name.all[i]=="nlp_log"){browser()}
        if(n.clust>1){
            tmpdiff = mean(avgcount[tmpfit@cluster==2]) - mean(avgcount[tmpfit@cluster==1])
            tmpind =  as.numeric( (tmpdiff > 0) + 1 )
            qprobtmp = qnorm(posterior(tmpfit))
            qprob = qprobtmp[,tmpind]
            qprob[is.infinite(qprob)] = -1*qprobtmp[is.infinite(qprob),3-tmpind]
            ### deal with some extreme clustering results ###
            if(sum(!is.infinite(qprob))>=2){
                qprob[qprob == Inf] = max(qprob[!is.infinite(qprob)])
                qprob[qprob == -Inf] = min(qprob[!is.infinite(qprob)])
            }else if(sum(!is.infinite(qprob))==1){
                if(qprob[!is.infinite(qprob)] >= 0){
                    qprob[qprob == Inf] = qprob[!is.infinite(qprob)]
                    qprob[qprob == -Inf] = qnorm(1/nrow(dat))
                }else{
                    qprob[qprob == Inf] = qnorm(1-1/nrow(dat))
                    qprob[qprob == -Inf] = qprob[!is.infinite(qprob)]
                }
            }else{
                qprob[qprob == Inf] = qnorm(1-1/nrow(dat))
                qprob[qprob == -Inf] = qnorm(1/nrow(dat))
            }
            #############

            qpost = data.frame("qprob" = qprob,
                               "group" = cbind(2-tmpfit@cluster,tmpfit@cluster-1)[,tmpind],
                               "prior" = tmpfit@prior[tmpind])

            qpost = as.matrix(qpost)

        }else{
            qpost = data.frame("qprob" = rep( qnorm(1-1/nrow(dat)), nrow(dat)),
                               "group" = rep(1,nrow(dat)),
                               "prior" = 1)
            qpost = as.matrix(qpost)
        }

        prior.all = c(prior.all,qpost[1,"prior"])
        post.all = cbind(post.all,qpost[,"qprob"])
        group.all = cbind(group.all, qpost[,"group"])
    }

    colnames(post.all) = name.all
    final.score = (apply(pnorm(post.all),1,mean) + pnorm(apply(post.all,1,mean)))/2

    cluster.ori = kmeans(cbind(dat[,vname.log],post.all[,vname.log]),centers=2) ### standardize?
    class.pos = 1 + 1*(cor(cluster.ori$cluster==2, avgcount)>0)
    prev.est0 = mean(cluster.ori$cluster==class.pos)
    prev.est = (mean(final.score)+prev.est0)/2

    final.score = Prob.Scale(final.score, prev.est)
    cut.MAP = quantile(final.score,1-prev.est)
    list("scores" = final.score, "cut.MAP" = cut.MAP)
}


#' Internal function to check parameters
#' @noRd
check_para <- function(dat.icd = NULL, dat.nlp = NULL, dat.note = NULL,
                       nm.phe = NULL, nm.ID = NULL, nm.utl = NULL,
                       p.icd = 0.001, n.icd = 100, p.nlp = 0.001, n.nlp = 100,
                       wgt.icd.nlp = c("con.5", "con1", "sd")[1],
                       yes.con = c(FALSE, TRUE)[1]){

    threshold = c(p.icd, n.icd, p.nlp, n.nlp)
    if((!is.numeric(threshold)) | (length(threshold)<4)){
        stop("Provide valid filter values to run MAP!")
    }

    if(is.null(nm.phe)){
        stop("Provide names of phenotypes (nm.phe) to run MAP!")
    }

    if(is.null(nm.ID)){
        stop("Provide variable name of ID (nm.ID) for merging the ICD, NLP, and note count data!")
    }

    if(is.null(nm.utl)){
        stop("Provide variable name of health utilizaitons (nm.utl) to adjust for!")
    }

    if(is.null(dat.icd) | is.null(dat.nlp) | is.null(dat.note)){
        stop("At least one of the ICD, NLP, or note count data are missing!")
    }

    if(any(!is.logical(yes.con))){
        stop("yes.con has to be TRUE or FALSE!")
    }

    if(any(!is.element(wgt.icd.nlp,c("con.5", "con1", "sd"))) ){
        stop("wgt.icd.nlp has to be one of {con.5, con1, sd}!")
    }

}

#' @title MAP algorithm
#' @description  Main function to perform MAP algorithm to calculate predicted probabilities of positive phenotype for each patient
#' based on NLP and ICD counts adjusted for healthcare utilization.
#' @param dat.icd Data containing ICD codes counts and must have
#' nm.ID indicating the patient IDs.
#' @param dat.nlp Data containing NLP counts and must have
#' nm.ID indicating the patient IDs.
#' @param dat.note Data containing health utilization count and must have
#' nm.ID indicating the patient IDs.
#' @param nm.phe A vector of phenotypes whose results are desired.
#' @param nm.ID Variable name for the patient IDs that are shared by dat.icd,
#' dat.nlp, and dat.note.
#' @param nm.utl Variable name for the health utilization count in dat.note.
#' @param p.icd Threshold in percentage to check if ICD data is good to run MAP.
#' @param n.icd Threshold in count to check if ICD data is good to run MAP
#' (ICD has to pass both p.icd and n.icd threshold to be included ).
#' @param p.nlp Threshold in percentage to check if NLP data is good to run MAP.
#' @param n.nlp Threshold in count to check if NLP data is good to run MAP.
#' (NLP has to pass both p.nlp and n.nlp threshold to be included ).
#' @param wgt.icd.nlp Weighting scheme to create a new variable using ICD and NLP.
#' @param yes.con A logical variable indicating if concomitant is desired.
#' @param yes.nlp A logical variable indicating if patients having zero ICD, but nonzero
#' NLP are set aside.
#' @return Returns a list with following objects:
#' \item{MAP}{A data.frame with each column being the probabilities for each phenotype.}
#' \item{cut.MAP}{A vector of cutoff values for the probabilities for each phenotype.}
#' @references Integrating natural language processing for high-throughput Multimodal
#' Automated Phenotyping (MAP) with Application to PheWAS. Katherine P. Liao, Jiehuan Sun,
#' Tianrun A. Cai, Nicholas Link, Chuan Hong, Jie Huang, Jennifer Huffman, Jessica Gronsbell, Yichi Zhang,
#' Yuk-Lam Ho, Jacqueline Honerlaw, Lauren Costa, Victor Castro, Vivian Gainer, Shawn Murphy,
#' Christopher J. O’Donnell, J. Michael Gaziano, Kelly Cho, Peter Szolovits, Isaac Kohane,
#' Sheng Yu, and Tianxi Cai with the VA Million Veteran Program
#' @examples
#' ## simulate data to test the algorithm
#' ## see the vignette for more examples
#' ncase = 500 ; ncontrol = 2000
#' dd = c(rep(0,ncontrol), rep(1,ncase))
#' note = rpois(ncase+ncontrol, 20)
#' ICD = rpois(ncase+ncontrol, dd*3 + 0.1 + 0.1*log(note))
#' NLP = rpois(ncase+ncontrol, dd*3 + 0.1 + 0.1*log(note))
#' dat.icd = data.frame(ID=1:length(ICD), disease = ICD)
#' dat.nlp = data.frame(ID=1:length(ICD), disease = NLP)
#' dat.note = data.frame(ID=1:length(ICD), note = note)
#' nm.phe = 'disease' ; nm.ID = 'ID' ; nm.utl = 'note';
#' res = MAP_PheWAS_main(dat.icd = dat.icd, dat.nlp = dat.nlp, dat.note = dat.note,
#' nm.phe = nm.phe, nm.ID = nm.ID, nm.utl = nm.utl)
#' boxplot(res$MAP~dd)

MAP_PheWAS_main <- function(dat.icd = NULL, dat.nlp = NULL, dat.note = NULL,
                            nm.phe = NULL, nm.ID = NULL, nm.utl = NULL,
                            p.icd = 0.001, n.icd = 100, p.nlp = 0.001, n.nlp = 100,
                            wgt.icd.nlp = c("con.5", "con1", "sd")[1],
                            yes.con = c(FALSE, TRUE)[1], yes.nlp = c(TRUE, FALSE)[1]){

    ### checking the data inputs
    check_para(dat.icd, dat.nlp, dat.note,  nm.phe, nm.ID, nm.utl,
               p.icd, n.icd, p.nlp, n.nlp, wgt.icd.nlp, yes.con)

    Noset = setdiff(nm.phe, union(colnames(dat.icd), colnames(dat.nlp)))
    if(length(Noset)>0){
        warning(cat("The following phenotypes do not exist in any data: \n", Noset,"\n") )
    }
    nm.phe = setdiff(nm.phe, Noset)

    ID = sort(intersect( intersect(dat.icd[,nm.ID], dat.nlp[,nm.ID]),
                         dat.note[,nm.ID]))

    dat.icd = dat.icd[match(ID, dat.icd[,nm.ID]),]
    dat.nlp = dat.nlp[match(ID, dat.nlp[,nm.ID]),]
    dat.note = dat.note[match(ID, dat.note[,nm.ID]),]
    note = log(dat.note[,nm.utl] + 1.0)

    MAPscore = NULL
    CUTmap = NULL
    for(pheno in nm.phe){
        cat(pheno,"\n")
        if(is.element(pheno, colnames(dat.icd))){
            icd = dat.icd[,pheno]
        }else{icd = rep(0, nrow(dat.icd))}

        if(is.element(pheno, colnames(dat.nlp))){
            nlp = dat.nlp[,pheno]
        }else{nlp = rep(0, nrow(dat.nlp))}

        dat = data.frame(icd=icd, nlp=nlp)
        if(any(is.na(dat))){
            stop(paste("No Missing values are allowed! Please check ", pheno, "!\n"))
        }

        filter.cc = apply(dat, 2, sum)
        filter.pp = apply(dat, 2, mean)

        icd.ind = NULL
        nlp.ind = NULL

        MAP.p00 = rep(0,nrow(dat))

        ### determine whether ICD and/or NLP satisfies the filtering criteria.
        if(filter.cc['icd'] > n.icd & filter.pp['icd'] > p.icd){

            if(filter.cc['nlp'] > n.nlp & filter.pp['nlp'] > p.nlp){
                if(yes.nlp){
                    icd.ind = (dat[,'icd'] > 0) | (dat[,'nlp'] > 0)
                }else{
                    icd.ind = (dat[,'icd'] > 0)
                }
                dat00 = dat[icd.ind,]
                note00 = note[icd.ind]

                ## choose weight to create a new variable based on ICD and NLP
                switch(wgt.icd.nlp,
                       "con.5"={tmpw = c(0.5,0.5)},
                       "con1"={tmpw = c(1,1)},
                       "sd"={tmpw = c(1/sd(dat00[,'icd']),1/sd(dat00[,'nlp']));
                       tmpw = tmpw/sum(tmpw)})
                dat00[,'icdnlp'] = round(as.matrix(dat00[,c('icd','nlp')])%*%tmpw)
                colnm = colnames(dat00)
                cat("**** Results w/ ICD+NLP ****\n")
            }else{
                icd.ind = (dat[,'icd'] > 0)
                dat00 = dat[icd.ind,]
                note00 = note[icd.ind]
                colnm = "icd"
                cat("Poor NLP Mapping for", pheno ,", Results w/ ICD only \n")
            }

        }else{

            if(filter.cc['nlp'] > n.nlp & filter.pp['nlp'] > p.nlp){
                nlp.ind = dat[,'nlp'] > 0
                dat00 = dat[nlp.ind,]
                note00 = note[nlp.ind]
                colnm = "nlp"
                cat("**** Results w/ NLP only ****\n")
            }else{
                colnm = NULL
                cat("Poor ICD and NLP Mapping for",pheno, ", use ICD >=1 as Classification \n")
            }
        }

        ### running MAP algorithm
        if(is.null(colnm)){
            MAP.p00 = 1.0*(dat[,'icd'] >= 1); cut.MAP = 0.5
        }else{
            dat00 = as.data.frame(cbind(dat00[,colnm], log(dat00[,colnm]+1), note=note00))
            colnames(dat00) = c(colnm, paste(colnm, "_log", sep=""), "note")
            junk = MAP_PheWAS_JS(dat = dat00, vname = colnm, yes.con = yes.con)
            cut.MAP = junk$cut.MAP
            if(is.null(icd.ind)){
                MAP.p00[nlp.ind] = junk$scores
            }else{
                MAP.p00[icd.ind] = junk$scores
            }
        }

        MAPscore = cbind(MAPscore, MAP.p00)
        CUTmap = c(CUTmap, cut.MAP)
    }

    rownames(MAPscore) = ID
    colnames(MAPscore) = nm.phe
    list(MAP = MAPscore, cut.MAP = CUTmap)
}

