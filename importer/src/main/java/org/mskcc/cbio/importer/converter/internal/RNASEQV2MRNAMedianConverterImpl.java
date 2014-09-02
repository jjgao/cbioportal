/** Copyright (c) 2014 Memorial Sloan-Kettering Cancer Center.
 *
 * This library is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY, WITHOUT EVEN THE IMPLIED WARRANTY OF
 * MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE.  The software and
 * documentation provided hereunder is on an "as is" basis, and
 * Memorial Sloan-Kettering Cancer Center 
 * has no obligations to provide maintenance, support,
 * updates, enhancements or modifications.  In no event shall
 * Memorial Sloan-Kettering Cancer Center
 * be liable to any party for direct, indirect, special,
 * incidental or consequential damages, including lost profits, arising
 * out of the use of this software and its documentation, even if
 * Memorial Sloan-Kettering Cancer Center 
 * has been advised of the possibility of such damage.
*/
package org.mskcc.cbio.importer.converter.internal;

import org.mskcc.cbio.importer.*;
import org.mskcc.cbio.importer.util.MapperUtil;
import org.mskcc.cbio.importer.model.*;

import org.apache.commons.logging.*;

import java.util.*;

/**
 * Class which implements the Converter interface for processing rna-seq (v2) RSEM files.
 */
public class RNASEQV2MRNAMedianConverterImpl extends RNASEQV2MRNAMedianConverterBase implements Converter {

	public RNASEQV2MRNAMedianConverterImpl(Config config, FileUtils fileUtils,
										   CaseIDs caseIDs, IDMapper idMapper)
	{
		this.config = config;
        this.fileUtils = fileUtils;
		this.caseIDs = caseIDs;
		this.idMapper = idMapper;
		this.LOG = LogFactory.getLog(RNASEQV2MRNAMedianConverterImpl.class);
		this.conversionType = ConversionType.TUMOR_ONLY;
	}

	public RNASEQV2MRNAMedianConverterImpl(Config config, FileUtils fileUtils,
	                                       CaseIDs caseIDs, IDMapper idMapper,
	                                       Log log, ConversionType conversionType)
	{
		this.config = config;
        this.fileUtils = fileUtils;
		this.caseIDs = caseIDs;
		this.idMapper = idMapper;
		this.LOG = log;
		this.conversionType = conversionType;
	}
}
