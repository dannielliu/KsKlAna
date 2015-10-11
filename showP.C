double pKs_eeKsK892(double ebeam)
{
  double mK0892 = 0.89594;
  double mKs = 0.497614;

  double EKs = 0.5*(ebeam - (pow(mK0892,2)-pow(mKs,2))/ebeam);
  double pKs = sqrt(pow(EKs,2)-pow(mKs,2));
  return pKs;
}

double pKs_eeKsKl(double ebeam)
{
  double mKs = 0.497614;
  double pKs = sqrt(pow(ebeam/2,2)-pow(mKs,2));
  return pKs;
}

int showP()
{
  double enes[21] = {2.0, 2.05, 2.10, 2.15, 2.175,
                     2.2, 2.2324, 2.3096, 2.3864, 2.396,
		     2.5, 2.6444, 2.6464, 2.7, 2.8,
		     2.9, 2.95, 2.981, 3.0, 3.02,
		     3.08};
  for (int i=0; i<21; i++){
    cout<<enes[i]<<"\t"<<pKs_eeKsKl(enes[i])<<"\t"<<pKs_eeKsK892(enes[i])<<"\t"<<pKs_eeKsKl(enes[i])-pKs_eeKsK892(enes[i])<<endl;
  }
  return 0;
}
