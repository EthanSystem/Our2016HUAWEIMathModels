package relief;

public class MainProcess
{
	public static void main(String args[])
	{
		String path="C:\\Users\\Administrator\\Desktop\\genotype.dat";
		String outpath="C:\\Users\\Administrator\\Desktop\\weight.txt";
		String pathid="C:\\Users\\Administrator\\Desktop\\WeiDianId.txt";
		String pathimportid="C:\\Users\\Administrator\\Desktop\\importantweidian.txt";
		Functions functions=new Functions();
		//functions.init(path);
		functions.CalMAF(0.05f, path,pathid);
		//functions.readImportweidian(pathimportid);
		//functions.mainProcedure(path,3, outpath, pathimportid);
		//functions.mainProcedure(path,1,outpath);
	}

}
