/** Copyright (c) 2012 Memorial Sloan-Kettering Cancer Center.
**
** This library is free software; you can redistribute it and/or modify it
** under the terms of the GNU Lesser General Public License as published
** by the Free Software Foundation; either version 2.1 of the License, or
** any later version.
**
** This library is distributed in the hope that it will be useful, but
** WITHOUT ANY WARRANTY, WITHOUT EVEN THE IMPLIED WARRANTY OF
** MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE.  The software and
** documentation provided hereunder is on an "as is" basis, and
** Memorial Sloan-Kettering Cancer Center 
** has no obligations to provide maintenance, support,
** updates, enhancements or modifications.  In no event shall
** Memorial Sloan-Kettering Cancer Center
** be liable to any party for direct, indirect, special,
** incidental or consequential damages, including lost profits, arising
** out of the use of this software and its documentation, even if
** Memorial Sloan-Kettering Cancer Center 
** has been advised of the possibility of such damage.  See
** the GNU Lesser General Public License for more details.
**
** You should have received a copy of the GNU Lesser General Public License
** along with this library; if not, write to the Free Software Foundation,
** Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA.
**/
package org.mskcc.cbio.portal.hotspots;

import java.util.Collection;

/**
 *
 * @author jgao
 */
public interface HotspotDetectiveParameters {
    
    /**
     * 
     * @return 
     */
    public Collection<Integer> getCancerStudyIds();
    
    /**
     * 
     * @param cancerStudyIds 
     */
    public void setCancerStudyIds(Collection<Integer> cancerStudyIds);
    
    /**
     * 
     * @return 
     */
    public Collection<String> getMutationTypes();
    
    /**
     * 
     * @param mutationTypes 
     */
    public void setMutationTypes(Collection<String> mutationTypes);

    /**
     * 
     * @return 
     */
    public Collection<Long> getEntrezGeneIds();
    
    /**
     * 
     * @param entrezGeneIds 
     */
    public void setEntrezGeneIds(Collection<Long> entrezGeneIds);

    /**
     * 
     * @return 
     */
    public Collection<Long> getExcludeEntrezGeneIds();
    /**
     * 
     * @param excludeEntrezGeneIds 
     */
    public void setExcludeEntrezGeneIds(Collection<Long> excludeEntrezGeneIds);

    /**
     * 
     * @return 
     */
    public int getThresholdHyperMutator();
    
    /**
     * 
     * @param thresholdHyperMutator 
     */
    public void setThresholdHyperMutator(int thresholdHyperMutator);
    
    /**
     * 
     * @return 
     */
    public int getThresholdSamples();
    
    /**
     * 
     * @param thresholdSamples 
     */
    public void setThresholdSamples(int thresholdSamples);
    
    /**
     * Use this for prefilter hotspots on single residue
     * @return 
     */
    public int getPrefilterThresholdSamplesOnSingleResidue();
    
    public void setPrefilterThresholdSamplesOnSingleResidue(int prefilterThresholdSamplesOnSingleResidue);
    
    /**
     * 
     * @return true if separate hotspots with different protein changes even on the same spot
     */
    public boolean getSeperateByProteinChangesForSingleResidueHotspot();
    
    /**
     * 
     * @param b 
     */
    public void setSeperateByProteinChangesForSingleResidueHotspot(boolean b);
    
    /**
     * 
     * @return 
     */
    public int getLinearSpotWindowSize();
    
    /**
     * 
     * @param linearSpotWindowSize 
     */
    public void setLinearSpotWindowSize(int linearSpotWindowSize);
    
    /**
     * 
     * @param  
     */
    public void setPtmHotspotWindowSize(int ptmHotspotWindowSize);
    
    /**
     * 
     * @param  
     */
    public int getPtmHotspotWindowSize();
    
    /**
     * 
     * @return 
     */
    public double getDistanceClosestAtomsThresholdFor3DHotspots();
    
    /**
     * 
     * @param anstrom 
     */
    public void setDistanceClosestAtomsThresholdFor3DHotspots(double anstrom);
    
    /**
     * 
     * @return 
     */
    public double getDistanceCAlphaThresholdFor3DHotspots();
    
    /**
     * 
     * @param anstrom 
     */
    public void setDistanceCAlphaThresholdFor3DHotspots(double anstrom);
    
    /**
     * 
     * @return 
     */
    public double getDistanceErrorThresholdFor3DHotspots();
    
    /**
     * 
     * @param anstrom 
     */
    public void setDistanceErrorThresholdFor3DHotspots(double anstrom);
    
    /**
     * 
     * @return 
     */
    public double getIdentpThresholdFor3DHotspots();
    
    /**
     * 
     * @param identpThreshold 
     */
    public void setIdentpThresholdFor3DHotspots(double identpThreshold);
    
    /**
     * 
     * @return 
     */
    public boolean getIncludingMismatchesFor3DHotspots();
    
    /**
     * 
     * @param includeMismatches 
     */
    public void setIncludingMismatchesFor3DHotspots(boolean includeMismatches);
    
    /**
     * 
     * @return 
     */
    public boolean getMergeOverlappingHotspots();
    
    /**
     * 
     * @param mergeOverlappingHotspots 
     */
    public void setMergeOverlappingHotspots(boolean mergeOverlappingHotspots);
}
