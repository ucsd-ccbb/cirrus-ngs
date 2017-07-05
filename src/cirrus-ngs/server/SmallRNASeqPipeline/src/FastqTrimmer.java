
public class FastqTrimmer {
	public void run(String readType, String shellFile, String inputFileName) {
		if (readType.equalsIgnoreCase("SE"))
		{
			FastqSETrimmer ft = new FastqSETrimmer();
			ft.setShellFile(shellFile);
				
			System.out.println("FastqTrimmer is processing:" + inputFileName);
			
			ft.m_fastqFile = inputFileName;
			
			ft.execute("", "");
		}
	}
}
