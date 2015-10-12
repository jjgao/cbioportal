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

import java.io.IOException;
import java.io.PrintWriter;
import java.util.HashMap;
import java.util.Map;
import javax.servlet.ServletException;
import javax.servlet.http.HttpServlet;
import javax.servlet.http.HttpServletRequest;
import javax.servlet.http.HttpServletResponse;
import org.apache.log4j.Logger;

/**
 *
 * @author jgao
 */
public class HotspotsServlet extends HttpServlet {
    private static Logger logger = Logger.getLogger(HotspotsServlet.class);
    public static final String HOTSPOT_TYPE = "hotspot_type";
    public static final String MUTATION_TYPE = "mutation_type";
    public static final String PTM_TYPE = "ptm_type";
    public static final String GENES = "genes";
    public static final String THRESHOLD_SAMPLES = "threshold_samples";
    public static final String THRESHOLD_MUTATIONS_HYPERMUTATOR = "threshold_hypermutator";
    public static final String PTM_HOTSPOT_WINDOW = "window_ptm";
    public static final String THRESHOLD_DISTANCE_CLOSEST_DISTANCE_CONTACT_MAP = "threshold_distance_closest_atom_3d";
    public static final String THRESHOLD_DISTANCE_C_ALPHA_CONTACT_MAP = "threshold_distance_c_alpha_3d";
    public static final String THRESHOLD_DISTANCE_ERROR_CONTACT_MAP = "threshold_distance_error_3d";
    public static final String THRESHOLD_UNIPROT_PDB_ALIGNMENT_IDENTP = "threshold_identp_3d";
    public static final String LINEAR_HOTSPOT_WINDOW = "window_linear";
    public static final String SINGLE_HOTSPOT_SEPARATE_BY_PROTEIN_CHANGE = "separate_by_protein_change";
    public static final String THRESHOLD_PREFILTER_RECURRENCE = "threshold_prefilter_recurrence";
    public static final String THRESHOLD_PVALUE = "threshold_pvalue";
    
    /** 
     * Processes requests for both HTTP <code>GET</code> and <code>POST</code> methods.
     * @param request servlet request
     * @param response servlet response
     * @throws ServletException if a servlet-specific error occurs
     * @throws IOException if an I/O error occurs
     */
    protected void processRequest(HttpServletRequest request, HttpServletResponse response)
            throws ServletException, IOException {
        String format = request.getParameter("format");
        
        if (format==null || format.equalsIgnoreCase("json")) {
            response.setContentType("application/json");
        }
        
        PrintWriter out = response.getWriter();
        try {
            HotspotMain.detectHotspot(getRequestParameterMap(request), out);
        } catch (HotspotException ex) {
            throw new ServletException(ex);
        } finally {            
            out.close();
        }
    }
    
    private Map<String, String> getRequestParameterMap(HttpServletRequest request) {
        Map<String, String[]> map = request.getParameterMap();
        Map<String, String> ret = new HashMap<String, String>();
        for (Map.Entry<String,String[]> entry : map.entrySet()) {
            ret.put(entry.getKey(), entry.getValue()[0]);
        }
        return ret;
    }

    public HotspotsServlet() {
    }
    
    // <editor-fold defaultstate="collapsed" desc="HttpServlet methods. Click on the + sign on the left to edit the code.">
    /** 
     * Handles the HTTP <code>GET</code> method.
     * @param request servlet request
     * @param response servlet response
     * @throws ServletException if a servlet-specific error occurs
     * @throws IOException if an I/O error occurs
     */
    @Override
    protected void doGet(HttpServletRequest request, HttpServletResponse response)
            throws ServletException, IOException {
        processRequest(request, response);
    }

    /** 
     * Handles the HTTP <code>POST</code> method.
     * @param request servlet request
     * @param response servlet response
     * @throws ServletException if a servlet-specific error occurs
     * @throws IOException if an I/O error occurs
     */
    @Override
    protected void doPost(HttpServletRequest request, HttpServletResponse response)
            throws ServletException, IOException {
        processRequest(request, response);
    }

    /** 
     * Returns a short description of the servlet.
     * @return a String containing servlet description
     */
    @Override
    public String getServletInfo() {
        return "Servlet to calculate and provide data of mutation hotspots";
    }// </editor-fold>
}
