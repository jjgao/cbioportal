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

import java.util.ArrayList;
import java.util.List;
import java.util.Set;
import java.util.TreeSet;
import org.apache.commons.lang.StringUtils;
import org.mskcc.cbio.portal.dao.DaoCancerStudy;
import org.mskcc.cbio.portal.dao.DaoGeneticProfile;
import org.mskcc.cbio.portal.model.CancerStudy;
import org.mskcc.cbio.portal.model.ExtendedMutation;

/**
 *
 * @author jgao
 */
public class HotspotImpl implements Hotspot {
    private MutatedProtein protein;
    private Set<Integer> residues;
    private List<ExtendedMutation> mutations;
    private List<Sample> samples;
    private String label;

    /**
     * 
     * @param gene
     * @param residues
     * @param label 
     */
    public HotspotImpl(MutatedProtein protein, Set<Integer> residues) {
        this.protein = protein;
        this.residues = residues;
        this.mutations = new ArrayList<ExtendedMutation>();
        this.samples = new ArrayList<Sample>();
    }
    
    /**
     * 
     * @return gene
     */
    @Override
    public MutatedProtein getProtein() {
        return protein;
    }
    
    /**
     * 
     * @return residues
     */
    @Override
    public Set<Integer> getResidues() {
        return residues;
    }
    
    @Override
    public List<ExtendedMutation> getMutations() {
        return mutations;
    }

    @Override
    public void addMutation(ExtendedMutation mutation) {
        mutations.add(mutation);
        CancerStudy cancerStudy = DaoCancerStudy.getCancerStudyByInternalId(
                DaoGeneticProfile.getGeneticProfileById(mutation.getGeneticProfileId()).getCancerStudyId());
        samples.add(new SampleImpl(mutation.getCaseId(), cancerStudy));
    }

    @Override
    public List<Sample> getSamples() {
        return samples;
    }

    /**
     * 
     * @param label 
     */
    public void setLabel(String label) {
        this.label = label;
    }
    
    /**
     * 
     * @return 
     */
    @Override
    public String getLabel() {
        if (label != null) {
            return label;
        }

        return protein.toString()+" "+StringUtils.join(new TreeSet<Integer>(getResidues()),";");
    }

    @Override
    public int hashCode() {
        return getLabel().hashCode();
    }

    @Override
    public boolean equals(Object obj) {
        if (obj == null) {
            return false;
        }
        if (!(obj instanceof Hotspot)) {
            return false;
        }
        final Hotspot other = (Hotspot) obj;
        return getLabel().equals(other.getLabel());
    }
}
