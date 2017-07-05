/**
 * The class is the entry to process alignment resulting files against genome or mirBase reference.
 * @author guorong
 *
 */
public class MiRNAParser {
	public static void main(String[] args) throws Exception {		
		/*
		if (args.length < 7)
		{
			System.out.println("##   Extract count matrix from SAM file aligned against hairpin genome.   ##");
			System.out.println("Usage: java -jar miRNAParser.jar counts_from_miRBase <Input_folder> <Output_folder> <Suffix> <gff_file> <design_table> <ignoreMismatchLocation>");
			System.out.println("");
			System.out.println("##   Extract count matrix from SAM file aligned against human genome.   ##");
			System.out.println("Usage: java -jar miRNAParser.jar counts_from_genome <Input_folder> <Output_folder> <Suffix> <gff_file> <design_table> <ignoreMismatchLocation>");
			return;
		}
		String packageName = args[0];
		String inputFile = args[1];
		String outputFile = args[2]; 
		String suffix = args[3]; 
		String gffFile = args[4];
		String designFile = args[5];
		String ignoreMismatchLocation = args[6];
*/

		String packageName = "counts_from_miRBase";
		String inputFolder = "/Users/guorongxu/Desktop/workspace/projects/jupyter-genomics_bitbucket/data/miRNASeq/";
		String outputFolder = "/Users/guorongxu/Desktop/workspace/projects/jupyter-genomics_bitbucket/data/miRNASeq/output/"; 
		String suffix = ".sam"; 
		String gffFile = "/Users/guorongxu/Desktop/workspace/projects/jupyter-genomics_bitbucket/data/miRNASeq/hsa.gff3";
		String designFile = "/Users/guorongxu/Desktop/workspace/projects/jupyter-genomics_bitbucket/data/miRNASeq/group.txt";
		String ignoreMismatchLocation = "true";

		MiRNAParser parser = new MiRNAParser();
		
		if ("counts_from_miRBase".equalsIgnoreCase(packageName))
			parser.parseMiRBase(inputFolder, outputFolder, suffix, gffFile, designFile, ignoreMismatchLocation);

		if ("counts_from_genome".equalsIgnoreCase(packageName))
			parser.parseGenome(inputFolder, outputFolder, suffix, gffFile, designFile, ignoreMismatchLocation);
	}
	
	public void parseMiRBase(String inputFolder, String outputFolder, String suffix, String gffFile, String designFile, 
			String ignoreMismatchLocation) throws Exception 
	{		
		System.out.println();
		System.out.println("System is parsing GFF: " + gffFile);
		System.out.println();
		
		GFFParser2 parser = new GFFParser2();
		parser.parseFile(gffFile);
		
		boolean ignore = Boolean.valueOf(ignoreMismatchLocation);
		
		MiRBaseCounter counter = new MiRBaseCounter();
		counter.setGFFParser(parser);
		counter.loadDesignTable(designFile);
		counter.setIgnoreMismatchLocation(ignore);
		counter.listAllFiles(inputFolder, suffix);

		System.out.println("No.\tFileName\tTotal\tUnmapped\tUndefined\t2Mismatches\t17Bases\t27Bases");
		
		for (int i = 0; i < counter.m_fileList.size(); i++)
		{
			counter.parseFile(counter.m_fileList.get(i), i);
			System.out.println((i + 1) + "\t" + counter.m_fileNameList.get(i) + "\t" + counter.m_totalReads + "\t" + counter.m_totalUnmappedReads 
					+ "\t" + counter.m_totalIgnoredReads + "\t" + counter.m_totalMismatch2 + "\t" + counter.m_totalIgnoredLess17 
					+ "\t" + counter.m_totalIgnoredLonger27);
		}
		
		MiRNAPrinter printer = new MiRNAPrinter(counter.m_fileNameList, counter.m_miRNACountsList, 
				counter.m_miRNAAnnotationList, counter.m_miRNAReadsList, counter.m_miRNAReadsIDList);
		
		if (!ignore)
		{
			printer.outputIsoformCountFile(outputFolder + "mirRNA.5p.isoform_counts.txt");
			printer.outputMiRNACountFile(outputFolder + "mirRNA.5p.statistics_counts.txt");
			printer.outputReadFile(outputFolder + "mirRNA.5p.isoforms.sam");
			printer.outputMiRNACountTable(outputFolder + "mirRNA.5p.isoform_counts_total.txt", true);
			printer.outputMiRNACountTable(outputFolder + "mirRNA.5p.isoform_counts_average.txt", false);
		}
		else
		{
			printer.outputIsoformCountFile(outputFolder + "mirRNA.all.isoform_counts.txt");
			printer.outputMiRNACountFile(outputFolder + "mirRNA.all.statistics_counts.txt");
			printer.outputReadFile(outputFolder + "mirRNA.all.isoforms.sam");
			printer.outputMiRNACountTable(outputFolder + "mirRNA.all.isoform_counts_total.txt", true);
			printer.outputMiRNACountTable(outputFolder + "mirRNA.all.isoform_counts_average.txt", false);
			
			MatrixSorter sorter = new MatrixSorter();
			sorter.sort(outputFolder + "mirRNA.all.isoform_counts_total.txt");
			sorter.sort(outputFolder + "mirRNA.all.isoform_counts_average.txt");
			
			CountPrinter cp = new CountPrinter();
			cp.execute(outputFolder);
			
			sorter.sort(outputFolder + "mirRNA.all.isoform_counts_total.normalized.txt");
			sorter.sort(outputFolder + "mirRNA.all.isoform_counts_average.normalized.txt");
			
			SAMPrinter sp = new SAMPrinter();
			sp.execute(outputFolder);
		}
	}
	
	public void parseGenome(String inputFile, String outputFile, String suffix, String gffFile, String designFile, 
			String ignoreMismatchLocation) throws Exception 
	{		
		GFFParser parser = new GFFParser();
		parser.parseFile(gffFile);
		
		boolean ignore = Boolean.valueOf(ignoreMismatchLocation);
		
		GenomeCounter counter = new GenomeCounter();
		counter.setGFFParser(parser);
		counter.loadDesignTable(designFile);
		counter.setIgnoreMismatchLocation(ignore);
		
		counter.listAllFiles(inputFile, suffix);
		
		System.out.println("No.\tFileName\tTotal\tUnmapped\tUndefined\t2Mismatches\t17Bases\t27Bases");
		
		for (int i = 0; i < counter.m_fileList.size(); i++)
		{
			counter.parseFile(counter.m_fileList.get(i), i);
			System.out.println((i+1) + "\t" + counter.m_fileNameList.get(i) + "\t" + counter.m_totalReads + "\t" + counter.m_totalUnmappedReads + "\t" + counter.m_totalIgnoredReads
					+ "\t" + counter.m_totalMismatch2 + "\t" + counter.m_totalIgnoredLess17 + "\t" + counter.m_totalIgnoredLonger27);
		}
		
		MiRNAPrinter printer = new MiRNAPrinter(counter.m_fileNameList, counter.m_miRNACountsList, 
				counter.m_miRNAAnnotationList, counter.m_miRNAReadsList, counter.m_miRNAReadsIDList);
		
		if (!ignore)
		{
			printer.outputIsoformCountFile(outputFile + "mirRNA.5p.isoform_counts.txt");
			printer.outputMiRNACountFile(outputFile + "mirRNA.5p.statistics_counts.txt");
			printer.outputReadFile(outputFile + "mirRNA.5p.isoforms.sam");
			printer.outputMiRNACountTable(outputFile + "mirRNA.5p.isoform_counts_total.txt", true);
			printer.outputMiRNACountTable(outputFile + "mirRNA.5p.isoform_counts_average.txt", false);
		}
		else
		{
			printer.outputIsoformCountFile(outputFile + "mirRNA.all.isoform_counts.txt");
			printer.outputMiRNACountFile(outputFile + "mirRNA.all.statistics_counts.txt");
			printer.outputReadFile(outputFile + "mirRNA.all.isoforms.sam");
			printer.outputMiRNACountTable(outputFile + "mirRNA.all.isoform_counts_total.txt", true);
			printer.outputMiRNACountTable(outputFile + "mirRNA.all.isoform_counts_average.txt", false);
		}
	}
}
