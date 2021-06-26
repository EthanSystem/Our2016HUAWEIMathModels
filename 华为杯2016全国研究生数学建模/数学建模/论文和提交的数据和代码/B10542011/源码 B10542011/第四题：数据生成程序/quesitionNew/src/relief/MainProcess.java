
package relief;

public class MainProcess
{
	public static void main(String args[])
	{
		String path="C:\\Users\\csp\\Desktop\\genotype.dat";
		//String outpath="C:\\Users\\Administrator\\Desktop\\weight.txt";
		//String pathid="C:\\Users\\Administrator\\Desktop\\WeiDianId.txt";
	//	String pathimportid="C:\\Users\\Administrator\\Desktop\\importantweidian.txt";
		Functions functions=new Functions();
		//functions.init(path);
		//functions.CalMAF(0.05f, path,pathid);
		//functions.readImportweidian(pathimportid);
		String qianzhui="C:\\Users\\csp\\Desktop\\question4\\orginal\\";
		String pathid="C:\\Users\\csp\\Desktop\\importantweidian\\";
		String outpath="C:\\Users\\csp\\Desktop\\result\\";
		int []xu={8};
		for(int i=0;i<1;i++)
		{
			String pathXZtype=qianzhui+xu[i]+".txt";
		    String newpathid=pathid+xu[i]+".txt";
		    String newoutpath=outpath+xu[i]+".txt";
			//functions.CalMAF(0.05f, path,newpathid,pathXZtype);
		   functions.mainProcedure(path,pathXZtype,1,newoutpath,newpathid);
		  //functions.mainProcedure(path,1,outpath);
		}
	}

}
