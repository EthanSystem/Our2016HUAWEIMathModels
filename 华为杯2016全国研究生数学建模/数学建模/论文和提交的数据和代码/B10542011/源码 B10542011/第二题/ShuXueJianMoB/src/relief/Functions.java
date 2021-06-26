package relief;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.*;
import tools.*;

public class Functions
{
	
	ArrayList<Sample> al=new ArrayList<Sample>();
	String []SNPname=new String[9445];//λ������
	//ArrayList<GeneMAF> MAF=new ArrayList<GeneMAF>();
	ArrayList<Integer> WeiDianId=new ArrayList<Integer>();
	int N=9445;//�������ص�
	int T=1000;//�����ĸ���/�������� ��һ���ѵ��������������ͬ��
	//��ά 
	public  void CalMAF(float yuzhi,String path,String pathid)
	{
		init(path);
		Sample head=al.get(0);
		String [] headName=head.getGenetype().split("\\s++");
		Sample tail=al.get(500);
		String [] tailName=head.getGenetype().split("\\s++");
		float MAFa;
		for(int i=0;i<N;i++)//λ��
		{
			
			int rate1=0;
			int rate2=0;
			GeneMAF gm=new GeneMAF();
			
			//ʶ���ʶ
			String gene=headName[i];
			char tag1=gene.charAt(0);
			 rate1++;
			 if(gene.charAt(1)==tag1)
				 rate1++;
			for(int j=1;j<500;j++)
			{
				Sample sample=al.get(j);
				 gene=sample.getGenetype().split("\\s++")[i];
				 if(gene.charAt(0)==tag1)
					 rate1++;
				 if(gene.charAt(1)==tag1)
					 rate1++;
			}
			String tai=tailName[i];
			char tag2=tai.charAt(0);
			    rate2++;
			if(gene.charAt(1)==tag2)
				 rate2++;
			for(int k=501;k<1000;k++)
			{
				
				 Sample sample=al.get(k);
				 gene=sample.getGenetype().split("\\s++")[i];
				 if(gene.charAt(0)==tag2)
					 rate2++;
				 if(gene.charAt(1)==tag2)
					 rate2++;
				
			}
			
			//ɸѡ����  
			if(tag1==tag2) //tag1��tag2��ʾλ����ͬ
			{
				
			   MAFa=Math.abs((rate1-rate2)/(float)1000);
			  if(MAFa>yuzhi)
			  {
				  WeiDianId.add(i);
				  System.out.println(MAFa);
			  }
			 
			}
			 //��ʾλ�㲻��ͬ
			else
			{
				 MAFa=Math.abs((1000-rate1-rate2)/(float)1000);
				 if(MAFa>yuzhi)
				  {
					  WeiDianId.add(i);
					  System.out.println(MAFa);
				  }
				
				
			}
			
			
			
		//	System.out.println("λ��:"+i);
		}
		
		
		writeWeidianId(pathid);
		
		System.out.println("ʣ��λ��:"+ WeiDianId.size());
	}
	  public void readImportweidian(String path)
	  {
		      FileInputStream fis=null;
			  InputStreamReader isr=null;
			  BufferedReader br=null;
			 try
			 {
				
				 fis=new FileInputStream(path);
				 isr=new InputStreamReader(fis);
				 br=new BufferedReader(isr);
				String a;
				
			     int id;
				 while((a=br.readLine())!=null)
				 {  
					 
					  id=Integer.parseInt(a);
					  WeiDianId.add(id);
					 
				 }
			
				 for(int i=0;i<WeiDianId.size();i++)
					 System.out.println(WeiDianId.get(i));
				 
			 }
			 catch(Exception e)
			 {
				 e.printStackTrace();
			 }
			 finally
			 {
				try {
					br.close();
				} catch (IOException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
				try {
					isr.close();
				} catch (IOException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
				try {
					fis.close();
				} catch (IOException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			 }
			 
	  }
	
	  public void writeWeidianId(String pathid)
	  {
		  
		  
		  FileWriter fw=null;
			 BufferedWriter bw=null;
			try {
				
				
				File a=new File(pathid);
				  if(!a.exists())
				  {
					  a.createNewFile();
				  }
				  fw=new FileWriter(a.getAbsoluteFile());
				  bw=new BufferedWriter(fw);
				
				  
				  for(int i=0;i<WeiDianId.size();i++)
				  {
					  bw.write(WeiDianId.get(i).toString());
					  bw.newLine();
				  }
				
			} 
			catch (Exception e) {
				// TODO: handle exception
			}
			
			finally
			{
				try {
					bw.close();
				} catch (IOException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
				try {
					fw.close();
				} catch (IOException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}
	  }
	//����
	Comparator<Neighbor> comparator = new Comparator<Neighbor>(){

		public int compare(Neighbor o1, Neighbor o2)
		{
			// TODO Auto-generated method stub
			if(o1.getWeight()<o2.getWeight())
        	{
        		return 1;
        	}
        	else if(o1.getWeight()==o2.getWeight())
        	{
        		return 0;
        	}
        	else
        		return -1;
		
		
		}};
	 //��ȡͬһ���Ľ���id
	public ArrayList<Neighbor> getSameNearstId(Sample training)
	{
		
		ArrayList<Neighbor> Neighbors=new ArrayList<Neighbor>();
		 int id=training.getId()-1;
		// System.out.println("id"+id);
		
		 if(training.isTag())//��ʾѵ�������в���״̬
		 {
			 Neighbors = lastCase(training, Neighbors, id);
			 //������ weight����Դ�С����
			 Collections.sort(Neighbors,comparator);
			// paixu(Neighbors);
		 }
		 else
		 {
			 
			 Neighbors = forecase(training, Neighbors, id);
			// System.out.println(Neighbors.get(1).getId()+":"+Neighbors.get(1).getWeight());
			 //������ weight����Դ�С����
			 Collections.sort(Neighbors,comparator);
			// paixu(Neighbors);
		    // System.out.println("Neighbors:"+Neighbors.size());
			 
		 }
		return Neighbors;
	}
	//ǰ500��û��������
	private ArrayList<Neighbor> forecase(Sample training, ArrayList<Neighbor> Neighbors,
			int id) {
		  int sameGene=0;//��ͬ�Ļ�����
		 String []trainingSplit=training.getGenetype().split("\\s+");
		for(int i=0;i<500;i++)
		 {
			if(i!=id)
			{
				sameGene=0;
				Neighbor neighbor=new Neighbor();
				Sample sample=al.get(i);
				String []sampleSplit=sample.getGenetype().split("\\s+");
				for(int j=0;j<9445;j++)
				{
					if(trainingSplit[j].equals(sampleSplit[j]))
					{
						sameGene++;
					}
					
				}
		
				//��id��Ȩֵ���뵽�ھ���
				neighbor.setId(i);
				neighbor.setWeight(sameGene/(float)N);
				Neighbors.add(neighbor);
				//System.out.println(Neighbors.size());
			}
		 }
		return Neighbors;
	}
	//��500���в�������
	private ArrayList<Neighbor> lastCase(Sample training, ArrayList<Neighbor> Neighbors,
			int id) {
		int sameGene=0;//��ͬ�Ļ�����
		String []trainingSplit=training.getGenetype().split("\\s+");
		//System.out.println("trainingSplit:"+trainingSplit.length);
		for(int i=500;i<1000;i++)
		 {
		   if(i!=id)	
		   {
			sameGene=0;
			Neighbor neighbor=new Neighbor();
			//��ȡ����
			Sample sample=al.get(i);
			String []sampleSplit=sample.getGenetype().split("\\s+");
		//	System.out.println("samplesplit:"+sampleSplit.length+"id"+i);
			for(int j=0;j<sampleSplit.length;j++)
			{
				if(trainingSplit[j].equals(sampleSplit[j]))
				{
					sameGene++;
				}
				
			}
			//��id��Ȩֵ���뵽�ھ���
			neighbor.setId(i);
			neighbor.setWeight(sameGene/(float)N);
			Neighbors.add(neighbor);
		   }
		 }
		return Neighbors;
	}
	//��ȡ�ڲ�ͬ����еĽ���id
	public ArrayList<Neighbor> getDistinctNearstId(Sample training)
	{
		 
		ArrayList<Neighbor> Neighbors=new ArrayList<Neighbor>();
		 int id=training.getId()-1;
		// int sameGene=0;//��ͬ�Ļ�����
		 if(!training.isTag())//��ʾѵ�������޲���״̬
		 {
			 Neighbors = lastCase(training, Neighbors, id);
			 //������ weight����Դ�С����
			 Collections.sort(Neighbors,comparator);
			// paixu(Neighbors);
		 }
		 //�в���״̬
		 else
		 {
			 
			 Neighbors = forecase(training, Neighbors, id);
			 //������ weight����Դ�С����
			 Collections.sort(Neighbors,comparator);
			// paixu(Neighbors);
			 
		 }
		return Neighbors;
	}
	
	
	
	/*private void paixu(ArrayList<Neighbor> Neighbors) {
		Collections.sort(Neighbors,new Comparator<Neighbor>() 
				 {
			        public int compare(Neighbor a,Neighbor b)
			        {
			        	if(a.getWeight()>b.getWeight())
			        	{
			        		return 1;
			        	}
			        	else if(a.getWeight()==b.getWeight())
			        	{
			        		return 0;
			        	}
			        	else
			        		return -1;
			        }
		
				 });
	}*/
	//д�ļ�
	
	public void writeTofile(String pathout,double weight[])
	{
		 FileWriter fw=null;
		 BufferedWriter bw=null;
		try {
			
			
			File a=new File(pathout);
			  if(!a.exists())
			  {
				  a.createNewFile();
			  }
			  fw=new FileWriter(a.getAbsoluteFile());
			  bw=new BufferedWriter(fw);
			
			  
			  for(int i=0;i<WeiDianId.size();i++)
			  {
				  bw.write(SNPname[WeiDianId.get(i)]+"\t"+weight[i]);
				  bw.newLine();
			  }
			
		} catch (Exception e) {
			// TODO: handle exception
		}
		
		finally
		{
			try {
				bw.close();
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			try {
				fw.close();
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
	}
	 
	//releif�㷨������,path��·����k�ǿ��Ʋ��ԵĽ��ڸ���
	public void mainProcedure(String path,int k,String outpath,String pathimportid)
	{
		init(path);
		readImportweidian(pathimportid);
		
	
		
		//��ʼ��һ��λ�������
		double weight[]=new double[WeiDianId.size()];
		float diffsame;
		float diffdistinct;
		for(int i=0;i<WeiDianId.size();i++)
		{
			weight[i]=0;
			
		}
		 int e=10;
		 int s=0;
		 int id;
		 int importantid;
		for(int j=0;j<1000;j++)//��������
		{
			
			  //���ѡȡһ���������з���
			 
				   
				  //id=(int)Math.round(Math.random()*(e-s)+s);
				 //  e=e+10;
				 //  s=s+10;
				 //  if(id==1000)
				//	   id=999;
			   
			  //int  id=(int) Math.round(Math.random()* 999);
			  //��ȡ����
			  Sample test=al.get(j);
			  ArrayList<Neighbor> sameNeighbors=getSameNearstId(test);
			  ArrayList<Neighbor> distinctNeighbors=getDistinctNearstId(test);
			  
			  for(int k1=0;k1<WeiDianId.size();k1++)
			  {
				  importantid=WeiDianId.get(k1);
				  for(int f=0;f<k;f++)
				  {
				  diffsame=diffValue(j, sameNeighbors.get(f).getId(),importantid);
				  diffdistinct=diffValue(j,distinctNeighbors.get(f).getId(),importantid);
				  weight[k1]=weight[k1]+(diffdistinct-diffsame)/(double)(N*T);
				  }
			  }
			  System.out.println(j);
		}
		
		writeTofile(outpath, weight);
		
		for(int t=0;t<10;t++)
		{
			System.out.println(SNPname[WeiDianId.get(t)]+":"+weight[t]);
		}
		
	}
  public int diffValue(int id1,int id2,int position)
  {
	   Sample a=al.get(id1);
	   Sample b=al.get(id2);
	  String []aGenetype=a.getGenetype().split("\\s+");
	  String []bGenetype=b.getGenetype().split("\\s+");
	  if(aGenetype[position].equals(bGenetype[position]))
	  {
		  return 0;
	  }
	  
	  
	  else
	  
		  return 1;
	  
	  
  }
	public void init(String path)
	{
	
		 FileInputStream fis=null;
		 InputStreamReader isr=null;
		 BufferedReader br=null;
		try
		{
			
			 fis=new FileInputStream(path);
			 isr=new InputStreamReader(fis);
			 br=new BufferedReader(isr);
			 String str="";
			 str=br.readLine();
			 SNPname=str.split("\\s++");
			 int id=1;
			 while((str=br.readLine())!=null)
			 {
				 Sample temp=new Sample();
				 if(id>500)
				 {
					 temp.setGenetype(str);
				     temp.setId(id);
				     temp.setTag(true);
				 }
				 else
				 {
					 temp.setGenetype(str);
				     temp.setId(id);
				     temp.setTag(false); 
				 }
				 id++;
				 
				 al.add(temp);
				 
			 }
			 //String[] training=al.get(250).getGenetype().split("\\s+");
			 //System.out.println(training.length);
			// System.out.println("����ٸ�:"+al.get(499).isTag());
			// System.out.println("����ٸ�:"+al.get(500).isTag());
			// System.out.println("����ٸ�:"+al.get(501).isTag());
			//System.out.println(al.size());
		}
		catch(Exception e)
		{
			e.printStackTrace();
		}
		
		
	}
	
	
	
}
	
	


