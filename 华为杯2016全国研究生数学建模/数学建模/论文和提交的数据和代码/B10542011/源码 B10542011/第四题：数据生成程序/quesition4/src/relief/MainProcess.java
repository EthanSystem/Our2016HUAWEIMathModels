
package relief;

public class MainProcess
{
	public static void main(String args[])
	{
		String path="C:\\Users\\csp\\Desktop\\genotype.dat";
		String outpath="C:\\Users\\csp\\Desktop\\weight.txt";
		//String pathid="C:\\Users\\Administrator\\Desktop\\WeiDianId.txt";
		String pathimportid="C:\\Users\\csp\\Desktop\\importantweidian.txt";
		Functions functions=new Functions();
		//functions.init(path);
		//functions.CalMAF(0.05f, path,pathid);
		//functions.readImportweidian(pathimportid);
		String qianzhui="C:\\Users\\csp\\Desktop\\question4\\orginal\\";
		String pathid="C:\\Users\\csp\\Desktop\\question4\\importantweidian\\";
		for(int i=10;i<=10;i++)
		{
			String pathXZtype=qianzhui+i+".txt";
			String newpathid=pathid+i+".txt";
			functions.CalMAF(0.05f, path,newpathid,pathXZtype);
		// functions.mainProcedure(path,pathXZtype,3, outpath, pathimportid);
		//functions.mainProcedure(path,1,outpath);
		}
	}

}
