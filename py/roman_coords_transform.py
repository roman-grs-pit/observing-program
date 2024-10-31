import json
import numpy as np
import matplotlib.patches as patches
from astropy.table import Table
import os
import matplotlib.pyplot as plt



##############################
# RomanCoordsTransform Class #
##############################
class RomanCoordsTransform:
    fpa_arrays: dict
    fpa_centers: dict
    fpa_zeros: dict
    wfi_arrays: dict
    wfi_centers: dict
    wfi_zeros: dict
    sca_arrays: dict
    sca_centers: dict
    sca_zeros: dict

    def __init__(
        self, file_path="./", json_file="roman_wfi_detector_coords.json"
    ):
        with open(os.path.join(file_path, json_file), "r") as read_file:
            data = read_file.read()
            json_dict = json.loads(data)

        self.fpa_arrays = json_dict["fpa_arrays"]
        self.fpa_centers = json_dict["fpa_centers"]
        self.fpa_zeros = json_dict["fpa_zeros"]

        self.wfi_arrays = json_dict["wfi_arrays"]
        self.wfi_centers = json_dict["wfi_centers"]
        self.wfi_zeros = json_dict["wfi_zeros"]

        self.sca_arrays = json_dict["sca_arrays"]
        self.sca_centers = json_dict["sca_centers"]
        self.sca_zeros = json_dict["sca_zeros"]

        # self.arrays = json_dict["wfi_local_arrays"]
        # self.array_centers = json_dict["wfi_local_centers"]
        # print(self.array_centers)

    def roman_fpa_layout(self, ax):
        self.roman_detector_layout(
            self.fpa_arrays, self.fpa_centers, ax, label="FPA", units="mm"
        )

    def roman_wfi_layout(self, ax):
        self.roman_detector_layout(
            self.wfi_arrays,
            self.wfi_centers,
            ax,
            label="WFI LOCAL",
            units="deg",
        )

    def roman_detector_layout(
        self, arrays, array_centers, ax, label="", units=""
    ):
        A_xmin = A_xmax = A_ymin = A_ymax = 0.0

        for i in range(18):
            sca = i + 1

            x0, y0 = array_centers["%i" % (sca)]

            A = arrays["%i" % (sca)]
            A_arr = np.array(A)

            Ax0 = np.min(A_arr[:, 0])
            Ax1 = np.max(A_arr[:, 0])
            Ay0 = np.min(A_arr[:, 1])
            Ay1 = np.max(A_arr[:, 1])

            if A_xmin > Ax0:
                A_xmin = Ax0
            if A_xmax < Ax1:
                A_xmax = Ax1
            if A_ymin > Ay0:
                A_ymin = Ay0
            if A_ymax < Ay1:
                A_ymax = Ay1

            ax.text(x0, y0, "%i" % sca)

            patch = patches.Polygon(A, ec="k", lw=1, fill=False)
            ax.add_patch(patch)

        Ax = A_xmax - A_xmin
        Ay = A_ymax - A_ymin
        dA = 0.05

        ax.scatter(0, 0, marker="x", s=100, c="k", lw=1.0)
        ax.scatter(0, 0, marker="o", s=200, fc="None", ec="k", lw=1.0)

        # print(A_xmin,A_xmax)
        # print(A_ymin,A_ymax)

        # print(A_xmin-dA*Ax,A_xmax+dA*Ax)
        # print(A_ymin-dA*Ay,A_ymax+dA*Ay)

        ax.set_aspect("equal")
        ax.set_xlim(A_xmin - dA * Ax, A_xmax + dA * Ax)
        ax.set_ylim(A_ymin - dA * Ay, A_ymax + dA * Ay)

        if units:
            units_str = "[%s]" % (units)
        else:
            units_str = ""

        xlabel_str = " ".join([label, "X", units_str])
        ylabel_str = " ".join([label, "Y", units_str])

        ax.set_xlabel(xlabel_str)
        ax.set_ylabel(ylabel_str)

    def inside_roman_fpa(self, x, y, ax):
        self.inside_roman_detector(self.fpa_arrays, x, y, ax)

    def inside_roman_wfi(self, x, y, ax):
        self.inside_roman_detector(self.wfi_arrays, x, y, ax)

    def inside_roman_detector(self, arrays, x, y, ax):
        for i in range(18):
            sca = i + 1

            # x0,y0 = array_centers[sca]

            A = arrays["%i" % (sca)]

            inside_filt = [
                RomanCoordsTransform.insidePolygon(A, [x0, y0])
                for x0, y0 in zip(x, y)
            ]

            x_in = x[inside_filt]
            y_in = y[inside_filt]

            ax.scatter(x_in, y_in, c="r", s=50)

    def measure_array_angles(self):
        for i in range(18):
            sca = i + 1
            print(sca)
            A = self.arrays["%i" % (sca)]
            ang = np.arctan2(A[0][1] - A[1][1], A[0][0] - A[1][0]) * 180 / np.pi
            print(ang)
            ang = np.arctan2(A[1][1] - A[2][1], A[1][0] - A[2][0]) * 180 / np.pi
            print(ang)
            ang = np.arctan2(A[2][1] - A[3][1], A[2][0] - A[3][0]) * 180 / np.pi
            print(ang)
            ang = np.arctan2(A[3][1] - A[0][1], A[3][0] - A[0][0]) * 180 / np.pi
            print(ang)

            # vec = np.arange(-360,360+90,90)
            # print( ang - vec)

            # if ang < -90 and > -180:
            # if ang < 0 and > -90:

            # if ang > 0 and < 90:
            # if ang > 90 and < 180:

            # print(array_centers[sca])
            # print()
            # x0,y0 = array_centers[sca]
            # for c1 in arrays[sca]:
            #    x1, y1 = c1
            #    ang = np.arctan2(y1 - y0, x1 - x0) * 180/np.pi
            #    print(ang)
            #    #ang = np.arctan2(y0 - y1, x0 - x1) * 180/np.pi
            #    #print(ang)
            # print()
            # for c1 in arrays[sca]:
            #    x1, y1 = c1
            #    #ang = np.arctan2(y1 - y0, x1 - x0) * 180/np.pi
            #    #print(ang)
            #    ang = np.arctan2(y0 - y1, x0 - x1) * 180/np.pi
            #    print(ang)
            print()
            print()

    def detector_to_fpa(self, xp, yp, sca):
        # print("Detector to FPA SCA %s" % (sca))

        """
        Transformation from detector [pixel] to FPA [mm] coordinates

        Parameters
        ----------
        xp : float
            Detector x-coordinates [pixel].

        yp : float
            Detector y-coordinates [pixel].

        sca : int or str
            Dectector number.

        Returns
        -------
        x : float
            FPA x-coordinates [mm].

        y : float
            FPA y-coordinates [mm].
        """

        # 1 mm = 1000 um

        # 10 um pixels
        # 10 um = 0.01 mm
        # 1 detector = 40.88 mm

        x0, y0 = self.fpa_zeros[str(sca)]
        # print(x0,y0)

        plate_scale_mm = 0.01  # mm/pixel

        xp *= plate_scale_mm
        yp *= plate_scale_mm

        x = xp + x0
        y = yp + y0

        return x, y

    def detector_to_wfi(self, xp, yp, sca):
        # print("Detector to WFI SCA %s" % (sca))
        """
        Transformation from detector [pixel] to WFI [deg] coordinates

        Parameters
        ----------
        xp : float
            Detector x-coordinates [pixel].

        yp : float
            Detector y-coordinates [pixel].

        sca : int or str
            Dectector number.

        Returns
        -------
        x : float
            WFI x-coordinates [deg].

        y : float
            WFI y-coordinates [deg].
        """

        x0, y0 = self.wfi_zeros[str(sca)]
        # print(x0,y0)

        plate_scale_arc = 0.11  # arcsec/pixel
        plate_scale_deg = plate_scale_arc / 3600.0  # deg/pixel

        xp *= plate_scale_deg
        yp *= -1 * plate_scale_deg

        x = xp + x0
        y = yp + y0

        return x, y

    def wfi_to_detector(self, x, y, sca):
        # print("WFI to detector SCA %s" % (sca))
        """
        Transformation from WFI [deg] to detector [pixel] coordinates

        Parameters
        ----------
        x : float
            WFI x-coordinates [deg].

        y : float
            WFI y-coordinates [deg].

        sca : int or str
            Dectector number.

        Returns
        -------
        xp : float
            Detector x-coordinates [pixel].

        yp : float
            Detector y-coordinates [pixel].
        """

        x0, y0 = self.wfi_zeros[str(sca)]
        # print(x0,y0)

        plate_scale_arc = 0.11  # arcsec/pixel
        plate_scale_deg = plate_scale_arc / 3600.0  # deg/pixel

        xp = x - x0
        yp = y - y0

        xp /= plate_scale_deg
        yp /= -1 * plate_scale_deg

        return xp, yp

    def fpa_to_detector(self, x, y, sca):
        # print("FPA to detector SCA %s" % (sca))
        """
        Transformation from FPA [mm] to detector [pixel] coordinates

        Parameters
        ----------
        x : float
            FPA x-coordinates [mm].

        y : float
            FPA y-coordinates [mm].

        sca : int or str
            Dectector number.

        Returns
        -------
        xp : float
            Detector x-coordinates [pixel].

        yp : float
            Detector y-coordinates [pixel].
        """
        # 1 mm = 1000 um

        # 10 um pixels
        # 10 um = 0.01 mm
        # 1 detector = 40.88 mm

        x0, y0 = self.fpa_zeros[str(sca)]
        # print(x0,y0)

        plate_scale_mm = 0.01  # mm/pixel

        xp = x - x0
        yp = y - y0

        xp /= plate_scale_mm
        yp /= plate_scale_mm

        return xp, yp

    def transform_detector_to_fpa(self, xp, yp, sca):
        """
        Transformation from detector [pixel] to FPA [mm] coordinates

        Parameters
        ----------
        xp : list or array
            Detector x-coordinates [pixel].

        yp : list or array
            Detector y-coordinates [pixel].

        sca : list
            Dectector number.

        Returns
        -------
        x : array
            FPA x-coordinates [mm].

        y : array
            FPA y-coordinates [mm].
        """

        if type(xp) == list:
            xp = np.array(xp)

        if type(yp) == list:
            yp = np.array(yp)

        x = []
        y = []
        for i in range(len(sca)):
            x0, y0 = self.detector_to_fpa(xp[i], yp[i], str(sca[i]))
            x.append(x0)
            y.append(y0)

        x = np.array(x)
        y = np.array(y)

        return x, y

    def transform_fpa_to_detector(self, x, y):
        """
        Transformation from FPA [mm] to detector [pixel] coordinates

        Parameters
        ----------
        x : list or array
            FPA x-coordinates [mm].

        y : list or array
            FPA y-coordinates [mm].

        Returns
        -------
        xp : array
            Detector x-coordinates [pixel].

        yp : array
            Detector y-coordinates [pixel].

        sca : list
            Dectector number.
        """

        if type(x) == list:
            x = np.array(x)

        if type(y) == list:
            y = np.array(y)

        xp = []
        yp = []
        sca = []

        for j in range(len(x)):
            for i in range(18):
                sca0 = i + 1

                A = self.fpa_arrays["%i" % (sca0)]

                if RomanCoordsTransform.insidePolygon(A, [x[j], y[j]]):
                    x0, y0 = self.fpa_to_detector(x[j], y[j], str(sca0))

                    xp.append(x0)
                    yp.append(y0)
                    sca.append(sca0)

        xp = np.array(xp)
        yp = np.array(yp)

        return xp, yp, sca

    def transform_detector_to_wfi(self, xp, yp, sca):
        """
        Transformation from detector [pixel] to WFI [deg] coordinates

        Parameters
        ----------
        xp : list or array
            Detector x-coordinates [pixel].

        yp : list or array
            Detector y-coordinates [pixel].

        sca : list
            Dectector number.

        Returns
        -------
        x : array
            WFI x-coordinates [deg].

        y : array
            WFI y-coordinates [deg].
        """

        if type(xp) == list:
            xp = np.array(xp)

        if type(yp) == list:
            yp = np.array(yp)

        x = []
        y = []
        for i in range(len(sca)):
            x0, y0 = self.detector_to_wfi(xp[i], yp[i], str(sca[i]))
            x.append(x0)
            y.append(y0)

        x = np.array(x)
        y = np.array(y)

        return x, y

    def transform_wfi_to_detector(self, x, y, clip_bad=0):
        """
        Transformation from WFI [deg] to detector [pixel] coordinates

        Parameters
        ----------
        x : list or array
            WFI x-coordinates [deg].

        y : list or array
            WFI y-coordinates [deg].

        Returns
        -------
        xp : array
            Detector x-coordinates [pixel].

        yp : array
            Detector y-coordinates [pixel].

        sca : list
            Dectector number.
        """

        if type(x) == list:
            x = np.array(x)

        if type(y) == list:
            y = np.array(y)

        xp = []
        yp = []
        sca = []

        if clip_bad:
            for j in range(len(x)):
                for i in range(18):
                    sca0 = i + 1

                    A = self.wfi_arrays["%i" % (sca0)]

                    if RomanCoordsTransform.insidePolygon(A, [x[j], y[j]]):
                        x0, y0 = self.wfi_to_detector(x[j], y[j], str(sca0))

                        xp.append(x0)
                        yp.append(y0)
                        sca.append(sca0)

        else:
            sep = np.zeros((18, x.shape[0]))

            for i in range(18):
                sca0 = i + 1
                # print("SCA =",sca)

                xc, yc = self.wfi_centers["%i" % (sca0)]

                sep[i, :] = np.sqrt((x - xc) ** 2 + (y - yc) ** 2)

            # print(sep)
            # print(np.min(sep,axis=0))
            ind = np.argmin(sep, axis=0)

            sca = [j + 1 for j in ind]

            for j in range(len(x)):
                x0, y0 = self.wfi_to_detector(x[j], y[j], str(sca[j]))

                xp.append(x0)
                yp.append(y0)

        xp = np.array(xp)
        yp = np.array(yp)

        return xp, yp, sca

    #############
    # prep data #
    #############
    @staticmethod
    def read_data(
        file_path="./", csv_file="WFIRST_WIM_190720.csv", col_prefix="col"
    ):
        tbl = Table.read(
            os.path.join(file_path, csv_file), data_start=17, format="ascii.csv"
        )
        print(f"column names:{tbl.colnames}")

        fields = [
            "SCA",
            "Wavelength",
            "Field",
            "SCA_X",
            "SCA_Y",
            "FPA_X",
            "FPA_Y",
            "WFI_LOCAL_X",
            "WFI_LOCAL_Y",
            "OPT_LOCAL_X",
            "OPT_LOCAL_Y",
            "OBS_FIELD_Y",
            "OBS_FIELD_Z",
        ]

        tbl.rename_column(col_prefix + "1", "SCA")
        tbl.rename_column(col_prefix + "2", "Wavelength")
        tbl.rename_column(col_prefix + "3", "Field")
        tbl.rename_column(col_prefix + "4", "SCA_X")
        tbl.rename_column(col_prefix + "5", "SCA_Y")
        tbl.rename_column(col_prefix + "6", "FPA_X")
        tbl.rename_column(col_prefix + "7", "FPA_Y")
        tbl.rename_column(col_prefix + "8", "WFI_LOCAL_X")
        tbl.rename_column(col_prefix + "9", "WFI_LOCAL_Y")
        tbl.rename_column(col_prefix + "10", "OPT_LOCAL_X")
        tbl.rename_column(col_prefix + "11", "OPT_LOCAL_Y")
        tbl.rename_column(col_prefix + "15", "OBS_FIELD_Y")
        tbl.rename_column(col_prefix + "16", "OBS_FIELD_Z")

        sub_tbl = tbl[fields]

        return sub_tbl

    @staticmethod
    def prep_data(
        file_path="./",
        csv_file="WFIRST_WIM_190720.csv",
        json_file="roman_wfi_detector_coords.json",
        col_prefix="col",
        write=0,
    ):
        tbl = RomanCoordsTransform.read_data(
            file_path=file_path, csv_file=csv_file, col_prefix=col_prefix
        )

        wfi_arrays = RomanCoordsTransform.build_roman_arrays(tbl)
        wfi_centers = RomanCoordsTransform.build_roman_array_centers(tbl)
        wfi_zeros = RomanCoordsTransform.build_roman_array_field(tbl, field=3)

        fpa_arrays = RomanCoordsTransform.build_roman_arrays(
            tbl, x_key="FPA_X", y_key="FPA_Y"
        )
        fpa_centers = RomanCoordsTransform.build_roman_array_centers(
            tbl, x_key="FPA_X", y_key="FPA_Y"
        )
        fpa_zeros = RomanCoordsTransform.build_roman_array_field(
            tbl, x_key="FPA_X", y_key="FPA_Y", field=3
        )

        sca_arrays = RomanCoordsTransform.build_roman_arrays(
            tbl, x_key="SCA_X", y_key="SCA_Y"
        )
        sca_centers = RomanCoordsTransform.build_roman_array_centers(
            tbl, x_key="SCA_X", y_key="SCA_Y"
        )
        sca_zeros = RomanCoordsTransform.build_roman_array_field(
            tbl, x_key="SCA_X", y_key="SCA_Y", field=3
        )

        json_dict = {}
        # json_dict["wfi_local_centers"] = array_centers
        # json_dict["wfi_local_arrays"] = arrays

        json_dict["fpa_arrays"] = fpa_arrays
        json_dict["fpa_centers"] = fpa_centers
        json_dict["fpa_zeros"] = fpa_zeros

        json_dict["wfi_arrays"] = wfi_arrays
        json_dict["wfi_centers"] = wfi_centers
        json_dict["wfi_zeros"] = wfi_zeros

        json_dict["sca_arrays"] = sca_arrays
        json_dict["sca_centers"] = sca_centers
        json_dict["sca_zeros"] = sca_zeros

        # print(json_dict)
        json_obj = json.dumps(json_dict)
        # print(json_obj)
        if write:
            print("Writing %s" % (json_file))
            with open(json_file, "w") as outfile:
                outfile.write(json_obj)

        return json_dict

    @staticmethod
    def build_roman_arrays(
        tbl, x_key="WFI_LOCAL_X", y_key="WFI_LOCAL_Y", wav=0.977
    ):
        arrays = {}
        for i in range(18):
            sca = i + 1

            filt_arr = (
                (tbl["SCA"] == sca)
                & (tbl["Field"] > 1)
                & (tbl["Wavelength"] == wav)
            )
            sub_tbl = tbl[filt_arr]

            A = np.array([sub_tbl[x_key], sub_tbl[y_key]]).T

            arrays["%i" % (sca)] = A.tolist()

        return arrays

    @staticmethod
    def build_roman_array_centers(
        tbl, x_key="WFI_LOCAL_X", y_key="WFI_LOCAL_Y", wav=0.977
    ):
        array_centers = {}
        for i in range(18):
            sca = i + 1

            filt = (
                (tbl["SCA"] == sca)
                & (tbl["Field"] == 1)
                & (tbl["Wavelength"] == wav)
            )
            x0, y0 = tbl[x_key, y_key][filt][0]

            p = np.array([x0, y0])

            array_centers["%i" % (sca)] = p.tolist()

        return array_centers

    @staticmethod
    def build_roman_array_field(
        tbl, x_key="WFI_LOCAL_X", y_key="WFI_LOCAL_Y", wav=0.977, field=1
    ):
        array_centers = {}
        for i in range(18):
            sca = i + 1

            filt = (
                (tbl["SCA"] == sca)
                & (tbl["Field"] == field)
                & (tbl["Wavelength"] == wav)
            )
            x0, y0 = tbl[x_key, y_key][filt][0]

            p = np.array([x0, y0])

            array_centers["%i" % (sca)] = p.tolist()

        return array_centers

    ##############################################

    # http://en.wikipedia.org/wiki/Point_in_polygon

    # Determining if a point lies on the interior of a polygon
    # Written by Paul Bourke
    # http://paulbourke.net/geometry/insidepoly/

    # Sum the angles of the point with the vertices of the polygon
    # Solution from Philippe Reverdy

    ############################################
    @staticmethod
    def insidePolygon(A, p):
        # format:

        #  points in polygon
        #  A = [[x0,y0],
        #       [x1,y1],
        #       [x2,y2]]

        #  test point
        #  p = [ x, y ]

        n = len(A)
        angle = 0
        for i in np.arange(n):
            x1 = A[i][0] - p[0]
            y1 = A[i][1] - p[1]
            x2 = A[(i + 1) % n][0] - p[0]
            y2 = A[(i + 1) % n][1] - p[1]
            angle += RomanCoordsTransform.angle2d(x1, y1, x2, y2)
        if np.abs(angle) < np.pi:
            return False
        else:
            return True

    @staticmethod
    def angle2d(x1, y1, x2, y2):
        theta1 = np.arctan2(y1, x1)
        theta2 = np.arctan2(y2, x2)
        dtheta = theta2 - theta1
        while dtheta > np.pi:
            dtheta -= 2 * np.pi
        while dtheta < -np.pi:
            dtheta += 2 * np.pi
        return dtheta

    def wfi_pa_rotation(self, PA):
    
        # PA = 0, theta=90
        # cartessian to PA -- PA to pointing to boresight
    
        #theta = (PA + 90.)
        theta = PA + 180

        #print(theta)
        theta *= np.pi/180. # radians
        #print(theta)

        R = np.array([
            [np.cos(theta), -1*np.sin(theta)],
            [np.sin(theta), np.cos(theta)],
        ])
    
        #R = np.array([
        #    [np.cos(theta), np.sin(theta)],
        #    [-1*np.sin(theta), np.cos(theta)],
        #])
    
        return R

    def wfi_sky_pointing(self, RA, Dec, PA, ax=None, alpha=1.0, c="k", ds9=True):

        if ax == None:
            fig = plt.figure()
            ax = fig.add_subplot(111)
    
        ax.scatter(RA, Dec, marker="o", s=50, fc="None", ec="tab:red")

        rot_arr = self.wfi_pa_rotation(PA)
        #print(rot_arr)
        #print()
        #print()
        if ds9:
            print("fk5")
            
            
        info = {'rot': rot_arr, 'PA': PA}
        
        for i in range(18):
        #for i in [0]:
            sca = i + 1
            A_vertices = np.array(self.wfi_arrays["%i" % (sca)])
            A_centers = np.array(self.wfi_centers["%i" % (sca)])
            #print(A)
            A_vertices[:,0] *= -1 # flip along y-axis (if SCAs are on the wrong side on sky)
            A_centers[0] *= -1 # flip along y-axis (if SCAs are on the wrong side on sky)
            #A_vertices[:,1] *= -1 # flip along x-axis (assume PA=0 is pointing toward the boresight)
            #A_centers[1] *= -1 # flip along x-axis (assume PA=0 is pointing toward the boresight)
        
            A_vert_rot = []
            for Av in A_vertices:
                Av_rot = np.dot(Av, rot_arr)   
                A_vert_rot.append(Av_rot)
                
            A_vert_rot = np.array(A_vert_rot)
            
            x0, y0 = np.dot(A_centers, rot_arr)             
            
            A_vert_rot[:,1] += Dec
            A_vert_rot[:,0] *= 1. / np.cos(A_vert_rot[:,1] * np.pi / 180.)
            A_vert_rot[:,0] += RA
            
            ax.text(x0/np.cos((y0+Dec)*np.pi/180.) + RA, y0 + Dec, "%i" % sca, c=c, alpha=alpha)
            patch = patches.Polygon(A_vert_rot, ec=c, lw=1, fill=False, alpha=alpha)
            ax.add_patch(patch)

            
            info[sca] = {
                "ra_cen": x0/np.cos((y0+Dec)*np.pi/180.) + RA, 
                "dec_cen": y0 + Dec, 
                "centers": [x0, y0], 
                "vertices": A_vert_rot,
            }
            
            if ds9:
                ###################
                # for ds9 regions #
                ###################
                a_list = []
                for a in A_vert_rot:
                    a_list += a.tolist() 
                #print(a_list)
            
                a_str = "polygon"
                for a in a_list:
                    a_str += " %.10f" % (a)
                a_str += " # text={%i}" % (sca)
                print(a_str)
                #print()
                ####################
            
        
        aspect = 1./np.cos(Dec * np.pi/180.)
        ax.set_aspect(aspect)
    
        x0, x1 = ax.get_xlim()
        if x1 > x0:
            ax.set_xlim(x1, x0)
        
        # in RA, Dec coordinates
        return info, ax

    def point_WFI_to_obj(self, RA, Dec, x_pix, y_pix, sca, ax=None, alpha=1.0, c="k", ds9=True):
        
        info = {}
        
        if ax == None:
            fig = plt.figure()
            ax = fig.add_subplot(111)
    
        #rct = RomanCoordsTransform(file_path=input_path)
        dx, dy = self.detector_to_wfi(x_pix, y_pix, sca) # WFI coords [deg]
        print(dx, dy)
    
        dec_ptg = Dec + dy
        ra_ptg = RA - dx/np.cos(dec_ptg*np.pi/180.)
        print(ra_ptg, dec_ptg)
    
        PA = 0. 
        roman_arrays, _ = self.wfi_sky_pointing(ra_ptg, dec_ptg, PA, ax=ax, alpha=alpha, c=c, ds9=ds9)
        info["roman_arrays"] = roman_arrays
    
        return ra_ptg, dec_ptg, ax, info
