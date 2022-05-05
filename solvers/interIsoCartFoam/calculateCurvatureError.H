const volScalarField K=interfaceForce.K();
const scalar Kref=2.5;
scalar err2=0.;
scalar errMax=0.;
scalar count=0.;
forAll(K,iCell)
{
  if (K[iCell]>0)
    {
      err2+=Foam::sqr((K[iCell]-Kref)/Kref);
      count+=1.;
      scalar err=Foam::mag(K[iCell]-Kref);
      if (err>errMax)
	errMax=err;
    }
}
if (count>0)
  err2/=count;
errMax/=Kref;
err2=Foam::sqrt(err2);
Info<<"err2="<<err2<<endl;
Info<<"errMax="<<errMax<<endl;
  
