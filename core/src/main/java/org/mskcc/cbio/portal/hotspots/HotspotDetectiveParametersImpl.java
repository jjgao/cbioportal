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
public class HotspotDetectiveParametersImpl implements HotspotDetectiveParameters {
    private Collection<Integer> cancerStudyIds;
    private Collection<String> mutationTypes; // missense or truncating
    private Collection<Long> entrezGeneIds;
    private Collection<Long>  excludeEntrezGeneIds;
    private int thresholdHyperMutator = -1;
    private int thresholdSamples;
    private boolean seperateByProteinChangesForSingleResidueHotspot;
    private int linearSpotWindowSize;
    private int ptmHotspotWindowSize;
    private double distanceThresholdClosestAtomsFor3DHotspots;
    private double distanceThresholdCAlphaFor3DHotspots;
    private double distanceErrorThresholdFor3DHotspots;
    private double identpThresholdFor3DHotspots;
    private boolean includingMismatchesFor3DHotspots;
    private int prefilterThresholdSamplesOnSingleResidue;
    private boolean mergeOverlappingHotspots;
    private double pValueThreshold;

    @Override
    public Collection<Integer> getCancerStudyIds() {
        return cancerStudyIds;
    }

    @Override
    public void setCancerStudyIds(Collection<Integer> cancerStudyIds) {
        this.cancerStudyIds = cancerStudyIds;
    }

    @Override
    public Collection<String> getMutationTypes() {
        return mutationTypes;
    }

    @Override
    public void setMutationTypes(Collection<String> mutationTypes) {
        this.mutationTypes = mutationTypes;
    }

    @Override
    public Collection<Long> getEntrezGeneIds() {
        return entrezGeneIds;
    }

    @Override
    public void setEntrezGeneIds(Collection<Long> entrezGeneIds) {
        this.entrezGeneIds = entrezGeneIds;
    }

    @Override
    public Collection<Long> getExcludeEntrezGeneIds() {
        return excludeEntrezGeneIds;
    }

    @Override
    public void setExcludeEntrezGeneIds(Collection<Long> excludeEntrezGeneIds) {
        this.excludeEntrezGeneIds = excludeEntrezGeneIds;
    }

    @Override
    public int getThresholdHyperMutator() {
        return thresholdHyperMutator;
    }

    @Override
    public void setThresholdHyperMutator(int thresholdHyperMutator) {
        this.thresholdHyperMutator = thresholdHyperMutator;
    }

    @Override
    public int getThresholdSamples() {
        return thresholdSamples;
    }

    @Override
    public void setThresholdSamples(int thresholdSamples) {
        this.thresholdSamples = thresholdSamples;
    }

    @Override
    public int getLinearSpotWindowSize() {
        return linearSpotWindowSize;
    }

    @Override
    public void setLinearSpotWindowSize(int linearSpotWindowSize) {
        this.linearSpotWindowSize = linearSpotWindowSize;
    }

    @Override
    public double getDistanceClosestAtomsThresholdFor3DHotspots() {
        return distanceThresholdClosestAtomsFor3DHotspots;
    }

    @Override
    public void setDistanceClosestAtomsThresholdFor3DHotspots(double anstrom) {
        this.distanceThresholdClosestAtomsFor3DHotspots = anstrom;
    }

    @Override
    public double getDistanceCAlphaThresholdFor3DHotspots() {
        return distanceThresholdCAlphaFor3DHotspots;
    }

    @Override
    public void setDistanceCAlphaThresholdFor3DHotspots(double anstrom) {
        this.distanceThresholdCAlphaFor3DHotspots = distanceThresholdCAlphaFor3DHotspots;
    }

    @Override
    public double getDistanceErrorThresholdFor3DHotspots() {
        return distanceErrorThresholdFor3DHotspots;
    }

    @Override
    public void setDistanceErrorThresholdFor3DHotspots(double distanceErrorThresholdFor3DHotspots) {
        this.distanceErrorThresholdFor3DHotspots = distanceErrorThresholdFor3DHotspots;
    }

    @Override
    public double getIdentpThresholdFor3DHotspots() {
        return identpThresholdFor3DHotspots;
    }

    @Override
    public void setIdentpThresholdFor3DHotspots(double identpThresholdFor3DHotspots) {
        this.identpThresholdFor3DHotspots = identpThresholdFor3DHotspots;
    }

    @Override
    public boolean getIncludingMismatchesFor3DHotspots() {
        return includingMismatchesFor3DHotspots;
    }

    @Override
    public void setIncludingMismatchesFor3DHotspots(boolean includeMismatches) {
        this.includingMismatchesFor3DHotspots = includeMismatches;
    }

    @Override
    public boolean getSeperateByProteinChangesForSingleResidueHotspot() {
        return seperateByProteinChangesForSingleResidueHotspot;
    }

    @Override
    public void setSeperateByProteinChangesForSingleResidueHotspot(boolean b) {
        this.seperateByProteinChangesForSingleResidueHotspot = b;
    }

    @Override
    public int getPrefilterThresholdSamplesOnSingleResidue() {
        return prefilterThresholdSamplesOnSingleResidue;
    }

    @Override
    public void setPrefilterThresholdSamplesOnSingleResidue(int prefilterThresholdSamplesOnSingleResidue) {
        this.prefilterThresholdSamplesOnSingleResidue = prefilterThresholdSamplesOnSingleResidue;
    }

    @Override
    public boolean getMergeOverlappingHotspots() {
        return mergeOverlappingHotspots;
    }

    @Override
    public void setMergeOverlappingHotspots(boolean mergeOverlappingHotspots) {
        this.mergeOverlappingHotspots = mergeOverlappingHotspots;
    }

    @Override
    public void setPtmHotspotWindowSize(int ptmHotspotWindowSize) {
        this.ptmHotspotWindowSize = ptmHotspotWindowSize;
    }

    @Override
    public int getPtmHotspotWindowSize() {
        return ptmHotspotWindowSize;
    }

    @Override
    public double getPValueThreshold() {
        return pValueThreshold;
    }

    @Override
    public void setPValueThreshold(double pValueThreshold) {
        this.pValueThreshold = pValueThreshold;
    }
}
