f0.loglike<-function(xj,nj,cj,dj)
{
	lbeta(xj+cj,nj-xj+dj)-lbeta(cj,dj)
}
fup.loglike<-function(xj,nj,cj,dj,pj)
{
	log(1/(1-pj))+log(1-pbeta(pj,xj+cj,nj-xj+dj))+lbeta(xj+cj,nj-xj+dj)-lbeta(cj,dj)
}
fdown.loglike<-function(xj,nj,cj,dj,pj)
{
	log(1/pj)+log(pbeta(pj,xj+cj,nj-xj+dj))+lbeta(xj+cj,nj-xj+dj)-lbeta(cj,dj)
}
iASeqmotiffit<-function(exprs,studyid,repid,refid,K,iter.max=100,tol=1e-3)
{
	set.seed(1)
	########################################################################
	#data preparation
	########################################################################
	sitenum<-nrow(exprs)
	studynum<-length(unique(studyid))
	experiment.unique<-unique(cbind(studyid,repid))
	colnames(experiment.unique)<-c("study","rep")
	xij<-matrix(0,nrow=sitenum,ncol=nrow(experiment.unique))
	nij<-matrix(0,nrow=sitenum,ncol=nrow(experiment.unique))	
	#extract reference allele frequencies
	flag=1
	for(d in 1:studynum)
	{
		studylab=unique(studyid)[d]
		repseq=which(experiment.unique[,1]==studylab)
		dn<-length(repseq)
		for(j in 1:dn)
		{
			replab=experiment.unique[repseq,2][j]
			sampleid=which(studyid==studylab & repid==replab & refid==0)
			xij[,flag]=exprs[,sampleid]
			flag=flag+1
		}
	}
	#extract total count for each allele in each sample
	flag=1
	for(d in 1:studynum)
	{
		studylab=unique(studyid)[d]
		repseq=which(experiment.unique[,1]==studylab)
		dn<-length(repseq)
		for(j in 1:dn)
		{
			replab=experiment.unique[repseq,2][j]
			sampleid=which(studyid==studylab & repid==replab)
			nij[,flag]=rowSums(exprs[,sampleid])
			flag=flag+1
		}
	}

	########################################################################
	#parameter preparation
	########################################################################
	pai=rep(0,K+1)
	p0j=rep(0,ncol(xij))
	p00j=rep(0,ncol(xij))
	#need fitting later
	c0j<-rep(1,ncol(xij))
	d0j<-rep(1,ncol(xij))

	cd<-rep(1,studynum)
	dd<-rep(1,studynum)

	#initialization
	pai<-c(0.8,rep(0.2/K,K))
	q0<-runif(K*studynum)
	q1<-runif(K*studynum)
	q2<-runif(K*studynum)
	qsum=(q0+q1+q2)
	qup=matrix(q1/qsum,K,studynum)
	qdown=matrix(q2/qsum,K,studynum)
	

	for(d in 1:studynum)
	{
		studylab=unique(studyid)[d]
		repseq=which(experiment.unique[,1]==studylab)
		dn<-length(repseq)
		tempx=rep(0,sitenum)
		tempn=rep(0,sitenum)
		for(j in 1:dn)
		{
			sampleid=repseq[j]
			tempx=xij[,sampleid]
			tempn=nij[,sampleid]
			prob=(tempx)/(tempn)	
			mu=mean(prob[which(tempn!=0)])
			v=var(prob[which(tempn!=0)])
			p00j[sampleid]=mean(prob[which(tempn!=0)])
		}
	}
	
	for(d in 1:studynum)
	{
		studylab=unique(studyid)[d]
		repseq=which(experiment.unique[,1]==studylab)
		dn<-length(repseq)
		tempx=rep(0,sitenum)
		tempn=rep(0,sitenum)
		for(j in 1:dn)
		{
			sampleid=repseq[j]
			tempx=xij[,sampleid]
			tempn=nij[,sampleid]
			addterm=p00j[sampleid]
			prob=(tempx+addterm*2)/(tempn+2)	
			mu=mean(prob[which(tempn!=0)])
			v=var(prob[which(tempn!=0)])
			c0j[sampleid]=mu*(mu*(1-mu)/v-1)
			d0j[sampleid]=(1-mu)*(mu*(1-mu)/v-1)
			p0j[sampleid]=mean(prob[which(tempn!=0)])
		}
	}
	c0j=abs(c0j)
	d0j=abs(d0j)
	loglike.old<--1e10
	########################################################################
	#model fitting
	########################################################################
	loglike0<-list()
	loglikeup<-list()
	loglikedown<-list()
	for(i.iter in 1:iter.max)
	{
		if(i.iter%%50==0)
		{
			print(paste("We have run ",i.iter," iterations for K=",K, sep=""))
		}
		err<-tol+1
		for(d in 1:studynum)
		{
			temploglike0=rep(0,sitenum)
			temploglikeup=rep(0,sitenum)
			temploglikedown=rep(0,sitenum)
			studylab=unique(studyid)[d]
			repseq=which(experiment.unique[,1]==studylab)
			dn<-length(repseq)
			for(j in 1:dn)
			{
				sampleid=repseq[j]
				temploglike0<-temploglike0+f0.loglike(xij[,sampleid],nij[,sampleid],c0j[sampleid],d0j[sampleid])
				temploglikeup<-temploglikeup+fup.loglike(xij[,sampleid],nij[,sampleid],cd[d],dd[d],p0j[sampleid])
				temploglikedown<-temploglikedown+fdown.loglike(xij[,sampleid],nij[,sampleid],cd[d],dd[d],p0j[sampleid])
			}
			loglike0[[d]]<-temploglike0
			loglikeup[[d]]<-temploglikeup
			loglikedown[[d]]<-temploglikedown
		}
		clustlike<-matrix(0,sitenum,K+1)
		condlikeup<-list()
		condlikedown<-list()
		for(d in 1:studynum)
		{
			condlikeup[[d]]<-matrix(0,sitenum,K+1)
			condlikedown[[d]]<-matrix(0,sitenum,K+1)
		}
		
		#for ai=0 case
		for(d in 1:studynum)
		{
			clustlike[,1]=clustlike[,1]+loglike0[[d]]
		}
		clustlike[,1]=clustlike[,1]+log(pai[1])
		#for ai!=0 case
		templike<-matrix(0,sitenum,3)
		for(k in 2:(K+1))
		{
			for(d in 1:studynum)
			{
				templike[,1]<-log(qup[k-1,d])+loglikeup[[d]]
				templike[,2]<-log(1-qup[k-1,d]-qdown[k-1,d])+loglike0[[d]]
				templike[,3]<-log(qdown[k-1,d])+loglikedown[[d]]
				tempmax<-pmax(templike[,1],templike[,2],templike[,3])
				for(z in 1:3)
				{
					templike[,z]<-exp(templike[,z]-tempmax)
				}
				tempsum<-templike[,1]+templike[,2]+templike[,3]
				clustlike[,k]<-clustlike[,k]+tempmax+log(tempsum)				
				condlikeup[[d]][,k]<-templike[,1]/tempsum
				condlikedown[[d]][,k]<-templike[,3]/tempsum
			}
			clustlike[,k]<-clustlike[,k]+log(pai[k])
		}

		tempmax<-apply(clustlike,1,max)
		for(k in 1:(K+1))	
		{
			clustlike[,k]<-exp(clustlike[,k]-tempmax)
		}	
		tempsum<-apply(clustlike,1,sum)

		#update motif occurrence rate
		for(k in 1:(K+1))
		{
			clustlike[,k]<-clustlike[,k]/tempsum
		}
		pai.new<-(apply(clustlike,2,sum)+1)/(sitenum+K+1)
	
		#update motifs
		qup.new<-matrix(0,K,studynum)
		qdown.new<-matrix(0,K,studynum)
		for(k in 1:K)
		{
			clustpsum<-sum(clustlike[,k+1])
			for(d in 1:studynum)
			{
				qup.new[k,d]<-(sum(clustlike[,k+1]*condlikeup[[d]][,k+1])+1)/(clustpsum+3)
				qdown.new[k,d]<-(sum(clustlike[,k+1]*condlikedown[[d]][,k+1])+1)/(clustpsum+3)
			}
		}

		loglike.new<-(sum(tempmax+log(tempsum))+sum(log(pai.new))
					+sum(log(qup.new)+log(qdown.new)+log(1-qup.new-qdown.new)))/sitenum

		## evaluate convergence
		err.pai<-max(abs(pai.new-pai)/pai)
		err.qup<-max(abs(qup.new-qup)/qup)
		err.qdown<-max(abs(qdown.new-qdown)/qdown)
		err<-max(err.pai, err.qup,err.qdown)

		pai<-pai.new
		qup<-qup.new
		qdown<-qdown.new
		loglike.old<-loglike.new

		if(err<tol)
		{
			break
		}
	}

	for(d in 1:studynum)
	{
		condlikeup[[d]]<-matrix(0,sitenum,K+1)
		condlikedown[[d]]<-matrix(0,sitenum,K+1)
	}

	#for ai=0 case
	for(d in 1:studynum)
	{
		clustlike[,1]=clustlike[,1]+loglike0[[d]]
	}
	clustlike[,1]=clustlike[,1]+log(pai[1])
	#for ai!=0 case
	templike<-matrix(0,sitenum,3)
	for(k in 2:(K+1))
	{
		for(d in 1:studynum)
		{
			templike[,1]<-log(qup[k-1,d])+loglikeup[[d]]
			templike[,2]<-log(1-qup[k-1,d]-qdown[k-1,d])+loglike0[[d]]
			templike[,3]<-log(qdown[k-1,d])+loglikedown[[d]]
			tempmax<-pmax(templike[,1],templike[,2],templike[,3])
			for(z in 1:3)
			{
				templike[,z]<-exp(templike[,z]-tempmax)
			}
			tempsum<-templike[,1]+templike[,2]+templike[,3]
			clustlike[,k]<-clustlike[,k]+tempmax+log(tempsum)				
			condlikeup[[d]][,k]<-templike[,1]/tempsum
			condlikedown[[d]][,k]<-templike[,3]/tempsum
		}
		clustlike[,k]<-clustlike[,k]+log(pai[k])
	}

	tempmax<-apply(clustlike,1,max)
	for(k in 1:(K+1))	
	{
		clustlike[,k]<-exp(clustlike[,k]-tempmax)
	}	
	tempsum<-apply(clustlike,1,sum)
	
	for(k in 1:(K+1))
	{
		clustlike[,k]<-clustlike[,k]/tempsum
	}

	p.post=1-(clustlike[,1])

	minus=sum(log(pai))+sum(log(qup)+log(qdown)+log(1-qup-qdown))
	
	result<-list(p.post=p.post, motif.prior=pai,motif.qup=qup,motif.qdown=qdown,clustlike=clustlike,
			c0j=c0j,d0j=d0j,loglike=loglike.new*sitenum-minus)
}
ASErawfit<-function(exprs,studyid,repid,refid)
{
	########################################################################
	#data preparation
	########################################################################
	sitenum<-nrow(exprs)
	studynum<-length(unique(studyid))
	experiment.unique<-unique(cbind(studyid,repid))
	colnames(experiment.unique)<-c("study","rep")
	xij<-matrix(0,nrow=sitenum,ncol=nrow(experiment.unique))
	nij<-matrix(0,nrow=sitenum,ncol=nrow(experiment.unique))	
	#extract reference allele frequencies
	flag=1
	for(d in 1:studynum)
	{
		studylab=unique(studyid)[d]
		repseq=which(experiment.unique[,1]==studylab)
		dn<-length(repseq)
		for(j in 1:dn)
		{
			replab=experiment.unique[repseq,2][j]
			sampleid=which(studyid==studylab & repid==replab & refid==0)
			xij[,flag]=exprs[,sampleid]
			flag=flag+1
		}
	}
	#extract total count for each allele in each sample
	flag=1
	for(d in 1:studynum)
	{
		studylab=unique(studyid)[d]
		repseq=which(experiment.unique[,1]==studylab)
		dn<-length(repseq)
		for(j in 1:dn)
		{
			replab=experiment.unique[repseq,2][j]
			sampleid=which(studyid==studylab & repid==replab)
			nij[,flag]=rowSums(exprs[,sampleid])
			flag=flag+1
		}
	}
	#result<-list(xij=xij,nij=nij,design=experiment.unique)

	########################################################################
	#parameter preparation
	########################################################################
	p00d=rep(0,ncol(xij))
	p0d=rep(0,studynum)
	p0dz=rep(0,studynum)
	z<-matrix(0,sitenum,studynum)
	b<-matrix(0,sitenum,studynum)
	B<-matrix(0,sitenum,studynum)
	c0d<-rep(1,studynum)
	d0d<-rep(1,studynum)

	xid<-matrix(0,nrow=sitenum,ncol=studynum)
	nid<-matrix(0,nrow=sitenum,ncol=studynum)	

	#pool counts
	for(d in 1:studynum)
	{
		studylab=unique(studyid)[d]
		repseq=which(experiment.unique[,1]==studylab)
		dn<-length(repseq)
		tempx=rep(0,sitenum)
		tempn=rep(0,sitenum)
		for(j in 1:dn)
		{
			sampleid=repseq[j]
			tempx=xij[,sampleid]+tempx
			tempn=nij[,sampleid]+tempn	
		}
		xid[,d]=tempx
		nid[,d]=tempn
		prob=(tempx)/(tempn)
		p00d[d]=mean(prob[which(tempn!=0)])
	}

	for(d in 1:studynum)
	{
		studylab=unique(studyid)[d]
		repseq=which(experiment.unique[,1]==studylab)
		dn<-length(repseq)
		tempx=rep(0,sitenum)
		tempn=rep(0,sitenum)
		for(j in 1:dn)
		{
			sampleid=repseq[j]
			tempx=xij[,sampleid]+tempx
			tempn=nij[,sampleid]+tempn	
		}
		#naive z statistic
		p0dz[d]=mean(tempx/tempn,na.rm=TRUE)
		z[,d]=abs(tempx/tempn-p0dz[d])

		#naive Bayes statistic
		addterm=p00d[d]
		prob=(tempx+2*addterm)/(tempn+2)
		mu=mean(prob[which(tempn!=0)])
		v=var(prob[which(tempn!=0)])
		c0d[d]=mu*(mu*(1-mu)/v-1)
		d0d[d]=(1-mu)*(mu*(1-mu)/v-1)
		p0d[d]=mean(prob[which(tempn!=0)])
		b[,d]=abs((tempx+2*addterm)/(tempn+2)-p0d[d])
		B[,d]=abs((tempx+c0d[d])/(tempx+c0d[d]+d0d[d])-p0d[d])
	}
	result<-list(z=z,b=b,B=B,c0d=c0d,d0d=d0d,p0d=p0d,p0dz=p0dz)
}


iASeqmotif<-function(exprs,studyid,repid,refid,K,iter.max=100,tol=1e-3)
{
	fitresult<-list()
	for(i in 1:length(K))
	{
		fitresult[[i]]<-iASeqmotiffit(exprs,studyid,repid,refid,K[i],iter.max=iter.max,tol=tol)
	}
	if(K[1]>1)
	{
		studynum=ncol(fitresult[[1]]$motif.qup)
	}
	else
	{
		studynum=length(fitresult[[1]]$motif.qup)
	}

	sitenum=nrow(exprs)
	bic<-rep(0,length(K))
	loglike<-rep(0,length(K))
	for(i in 1:length(K))
	{
		loglike[i]<-fitresult[[i]]$loglike
	}
	for(i in 1:length(K))
	{
		bic[i]<--2*fitresult[[i]]$loglike+(K[i]+K[i]*studynum*2+2*ncol(exprs)/2)*log(sitenum)
	}
	bestflag=which(bic==min(bic))
	result<-list(bestmotif=fitresult[[bestflag]],bic=cbind(K,bic),loglike=cbind(K,loglike))
}

###########################################################
#Motif plot
###########################################################
plotMotif<-function(bestmotif,title="",cutoff)
{
	  par(mai=c(0.7,1,0.7,0.2))
	  layout(matrix(1:4,ncol=4))
          u<-1:dim(bestmotif$motif.qup)[2]
          v<-1:dim(bestmotif$motif.qup)[1]
          image(u,v,t(bestmotif$motif.qup),
          col=gray(seq(from=1,to=0,by=-0.1)),xlab="StudyID",yaxt = "n",
		ylab="Motifs",main=paste(title,"V",sep=""),
		cex.axis=2,cex.lab=2,cex.main=2)
	  axis(2,at=1:length(v),,cex.axis=2,cex.lab=3)
          for(i in 1:(length(u)+1))
          {
                abline(v=(i-0.5))
          }
          for(i in 1:(length(v)+1)) 
          {
                abline(h=(i-0.5))
          }
	  image(u,v,t(bestmotif$motif.qdown),
          col=gray(seq(from=1,to=0,by=-0.1)),xlab="StudyID",yaxt = "n",
		ylab="Motifs",main=paste(title,"W",sep=""),
		cex.axis=2,cex.lab=2,cex.main=2)
	  axis(2,at=1:length(v),cex.axis=2,cex.lab=3)
          for(i in 1:(length(u)+1))
          {
                abline(v=(i-0.5))
          }
          for(i in 1:(length(v)+1)) 
          {
                abline(h=(i-0.5))
          }
	  Ng=nrow(bestmotif$clustlike)
	  genecount=floor(bestmotif$motif.p[-1]*Ng)
	  NK=nrow(bestmotif$motif.qup)
	  plot(0,0.7,pch=".",xlim=c(0,1.2),ylim=c(0.75,NK+0.25),frame.plot=FALSE,axes=FALSE,xlab="No. of SNPs",ylab=""
		,cex.axis=1.5,cex.lab=2,cex.main=3,main=expression(pi))
	  segments(0,0.7,bestmotif$motif.p[2],0.7)
	  rect(0,1:NK-0.3,bestmotif$motif.p[-1],1:NK+0.3,col="dark grey")
	  mtext(1:NK,at=1:NK,side=2,cex=1.5)
	  text(bestmotif$motif.p[-1]+0.55,1:NK,labels=floor(bestmotif$motif.p[-1]*Ng),cex=1)

	  clustnum=sapply(1:NK, function(k) 
		length(which(sapply(1:nrow(bestmotif$clustlike), function(i) which.max(bestmotif$clustlike[i,]))==(k+1) & 
									bestmotif$clustlike[,(k+1)]>cutoff)))
	  plot(0,0.7,pch=".",xlim=c(0,1.2),ylim=c(0.75,NK+0.25),frame.plot=FALSE,axes=FALSE,xlab="No. of SNPs",ylab=""
		,cex.axis=1.5,cex.lab=2,cex.main=3,main=expression(a[i]))
	  segments(0,0.7,clustnum[1]/Ng,0.7)
	  rect(0,1:NK-0.3,clustnum/Ng,1:NK+0.3,col="dark grey")
	  mtext(1:NK,at=1:NK,side=2,cex=1.5)
	  text(clustnum/Ng+0.55,1:NK,labels=clustnum,cex=1)
}

###########################################################
#BIC plot
###########################################################
plotBIC<-function(fitted_cormotif)
{
	plot(fitted_cormotif$bic[,1], fitted_cormotif$bic[,2], type="b", xlab="Non-null Motif Number", ylab="BIC", main="BIC")	
	lines(fitted_cormotif$bic[,1], fitted_cormotif$bic[,2])
}

singleEMfit<-function(exprs,studyid,repid,refid,iter.max=100,tol=1e-3)
{
	set.seed(1)
	########################################################################
	#data preparation
	########################################################################
	sitenum<-nrow(exprs)
	studynum<-length(unique(studyid))
	experiment.unique<-unique(cbind(studyid,repid))
	colnames(experiment.unique)<-c("study","rep")
	xij<-matrix(0,nrow=sitenum,ncol=nrow(experiment.unique))
	nij<-matrix(0,nrow=sitenum,ncol=nrow(experiment.unique))	
	#extract reference allele frequencies
	flag=1
	for(d in 1:studynum)
	{
		studylab=unique(studyid)[d]
		repseq=which(experiment.unique[,1]==studylab)
		dn<-length(repseq)
		for(j in 1:dn)
		{
			replab=experiment.unique[repseq,2][j]
			sampleid=which(studyid==studylab & repid==replab & refid==0)
			xij[,flag]=exprs[,sampleid]
			flag=flag+1
		}
	}
	#extract total count for each allele in each sample
	flag=1
	for(d in 1:studynum)
	{
		studylab=unique(studyid)[d]
		repseq=which(experiment.unique[,1]==studylab)
		dn<-length(repseq)
		for(j in 1:dn)
		{
			replab=experiment.unique[repseq,2][j]
			sampleid=which(studyid==studylab & repid==replab)
			nij[,flag]=rowSums(exprs[,sampleid])
			flag=flag+1
		}
	}
	#result<-list(xij=xij,nij=nij,design=experiment.unique)

	########################################################################
	#parameter preparation
	########################################################################
	p0j=rep(0,ncol(xij))
	p00j=rep(0,ncol(xij))
	#need fitting later
	c0j<-rep(1,ncol(xij))
	d0j<-rep(1,ncol(xij))

	cd<-rep(1,studynum)
	dd<-rep(1,studynum)

	#initialization
	qup=rep(0.1,studynum)
	qdown=rep(0.1,studynum)

	for(d in 1:studynum)
	{
		studylab=unique(studyid)[d]
		repseq=which(experiment.unique[,1]==studylab)
		dn<-length(repseq)
		tempx=rep(0,sitenum)
		tempn=rep(0,sitenum)
		for(j in 1:dn)
		{
			sampleid=repseq[j]
			tempx=xij[,sampleid]
			tempn=nij[,sampleid]
			prob=(tempx)/(tempn)	
			mu=mean(prob[which(tempn!=0)])
			v=var(prob[which(tempn!=0)])
			p00j[sampleid]=mean(prob[which(tempn!=0)])
		}
	}
	
	for(d in 1:studynum)
	{
		studylab=unique(studyid)[d]
		repseq=which(experiment.unique[,1]==studylab)
		dn<-length(repseq)
		tempx=rep(0,sitenum)
		tempn=rep(0,sitenum)
		for(j in 1:dn)
		{
			sampleid=repseq[j]
			tempx=xij[,sampleid]
			tempn=nij[,sampleid]
			addterm=p00j[sampleid]
			prob=(tempx+addterm*2)/(tempn+2)	
			mu=mean(prob[which(tempn!=0)])
			v=var(prob[which(tempn!=0)])
			c0j[sampleid]=mu*(mu*(1-mu)/v-1)
			d0j[sampleid]=(1-mu)*(mu*(1-mu)/v-1)
			p0j[sampleid]=mean(prob[which(tempn!=0)])
		}
	}
	c0j=abs(c0j)
	d0j=abs(d0j)
	
	loglike.old<--1e10
	########################################################################
	#model fitting
	########################################################################
	loglike0<-list()
	loglikeup<-list()
	loglikedown<-list()
	for(i.iter in 1:iter.max)
	{
		err<-tol+1
		for(d in 1:studynum)
		{
			temploglike0=rep(0,sitenum)
			temploglikeup=rep(0,sitenum)
			temploglikedown=rep(0,sitenum)
			studylab=unique(studyid)[d]
			repseq=which(experiment.unique[,1]==studylab)
			dn<-length(repseq)
			for(j in 1:dn)
			{
				sampleid=repseq[j]
				temploglike0<-temploglike0+f0.loglike(xij[,sampleid],nij[,sampleid],c0j[sampleid],d0j[sampleid])
				temploglikeup<-temploglikeup+fup.loglike(xij[,sampleid],nij[,sampleid],cd[d],dd[d],p0j[sampleid])
				temploglikedown<-temploglikedown+fdown.loglike(xij[,sampleid],nij[,sampleid],cd[d],dd[d],p0j[sampleid])
			}
			loglike0[[d]]<-temploglike0
			loglikeup[[d]]<-temploglikeup
			loglikedown[[d]]<-temploglikedown
		}
	
		#compute posterior cluster membership
		condlike<-list()
		for(d in 1:studynum)
		{
			condlike[[d]]<-matrix(0,sitenum,3)
		}
		
		templike<-matrix(0,sitenum,3)
		for(d in 1:studynum)
		{
			templike[,1]<-log(qup[d])+loglikeup[[d]]
			templike[,2]<-log(1-qup[d]-qdown[d])+loglike0[[d]]
			templike[,3]<-log(qdown[d])+loglikedown[[d]]
			tempmax<-pmax(templike[,1],templike[,2],templike[,3])
			for(z in 1:3)
			{
				templike[,z]<-exp(templike[,z]-tempmax)
			}
			tempsum<-templike[,1]+templike[,2]+templike[,3]				
			condlike[[d]][,1]<-templike[,1]/tempsum
			condlike[[d]][,2]<-templike[,2]/tempsum
			condlike[[d]][,3]<-templike[,3]/tempsum
		}
		
		qup.new<-rep(0,studynum)
		qdown.new<-rep(0,studynum)
		for(d in 1:studynum)
		{
			qup.new[d]<-(sum(condlike[[d]][,1])+1)/(sitenum+3)
			qdown.new[d]<-(sum(condlike[[d]][,3])+1)/(sitenum+3)
		}

		
		loglike.new<-(sum(tempmax+log(tempsum))
					+sum(log(qup.new)+log(qdown.new)+log(1-qup.new-qdown.new)))/sitenum
		
		## evaluate convergence
		err.qup<-max(abs(qup.new-qup)/qup)
		err.qdown<-max(abs(qdown.new-qdown)/qdown)
		err<-max(err.qup,err.qdown)

		qup<-qup.new
		qdown<-qdown.new
		loglike.old<-loglike.new

		if(err<tol)
		{
			break
		}
	}

	for(d in 1:studynum)
	{
		condlike[[d]]<-matrix(0,sitenum,3)
	}
	templike<-matrix(0,sitenum,3)
	for(d in 1:studynum)
	{
		templike[,1]<-log(qup[d])+loglikeup[[d]]
		templike[,2]<-log(1-qup[d]-qdown[d])+loglike0[[d]]
		templike[,3]<-log(qdown[d])+loglikedown[[d]]
		tempmax<-pmax(templike[,1],templike[,2],templike[,3])
		for(z in 1:3)
		{
			templike[,z]<-exp(templike[,z]-tempmax)
		}

		tempsum<-templike[,1]+templike[,2]+templike[,3]				
		condlike[[d]][,1]<-templike[,1]/tempsum
		condlike[[d]][,2]<-templike[,2]/tempsum
		condlike[[d]][,3]<-templike[,3]/tempsum

	}	
	p.post2=matrix(0,sitenum,studynum)
	for(d in 1:studynum)
	{
		p.post2[,d]=condlike[[d]][,2]
	}
	result<-list(p.study=p.post2,motif.qup=qup,motif.qdown=qdown,condlike=condlike)
}