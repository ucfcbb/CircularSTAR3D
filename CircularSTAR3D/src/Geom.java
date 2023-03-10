import java.util.*;

import org.ejml.data.*;
import org.ejml.interfaces.decomposition.SingularValueDecomposition;
import org.ejml.ops.*;
import org.ejml.factory.*;

public class Geom {
	public static Point centroid(List<Point> coord){
		double x=0., y=0., z=0.;
		for(Point P: coord){
			x+=P.x;
			y+=P.y;
			z+=P.z;
		}
		return new Point(x/coord.size(), y/coord.size(), z/coord.size());
	}
	
	public static Point centroid(DenseMatrix64F X){
		double x=0., y=0., z=0.;
		for(int i=0; i<X.numRows; i++){
			x+=X.get(i, 0);
			y+=X.get(i, 1);
			z+=X.get(i, 2);
		}
		return new Point(x/X.numRows, y/X.numRows, z/X.numRows);
	}
	
	public static double distance(double x1, double y1, double z1, double x2, double y2, double z2){
		return Math.sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2));
	}

	public static double RMSD(DenseMatrix64F X, DenseMatrix64F Y){
		if(X.numRows != Y.numRows) throw new RuntimeException("From RMSD: the lengths of X and Y are not same");
		
		double rmsd=0.;
		for(int i=0; i<X.numRows; i++){
			rmsd+=(X.get(i, 0)-Y.get(i, 0))*(X.get(i, 0)-Y.get(i, 0));
			rmsd+=(X.get(i, 1)-Y.get(i, 1))*(X.get(i, 1)-Y.get(i, 1));
			rmsd+=(X.get(i, 2)-Y.get(i, 2))*(X.get(i, 2)-Y.get(i, 2));
		}
		return Math.sqrt(rmsd/X.numRows);
	}
	
	public static DenseMatrix64F Kabsch(DenseMatrix64F X, DenseMatrix64F Y){
		DenseMatrix64F XT=new DenseMatrix64F(X.numCols, X.numRows);
		CommonOps.transpose(X, XT);
		DenseMatrix64F Z=new DenseMatrix64F(3,3);
		CommonOps.mult(XT, Y, Z);
		
		SingularValueDecomposition<DenseMatrix64F> svd = DecompositionFactory.svd(Z.numRows,Z.numCols,true,true,false);

		if(!svd.decompose(Z))
			throw new RuntimeException("Decomposition failed");
		DenseMatrix64F U = svd.getU(null,false);
        //DenseMatrix64F W = svd.getW(null);
        DenseMatrix64F V = svd.getV(null,false);
        
        double chi=Math.signum(CommonOps.det(Z));
        DenseMatrix64F C=new DenseMatrix64F(3,3);
        C.set(0,0,1.);
        C.set(1,1,1.);
        C.set(2,2,chi);
        
        CommonOps.transpose(U);
        
        DenseMatrix64F T=new DenseMatrix64F(3,3);
        DenseMatrix64F R=new DenseMatrix64F(3,3);
        CommonOps.mult(V, C, T);
        CommonOps.mult(T, U, R);
        return R;
	}

	public static DenseMatrix64F translation(DenseMatrix64F X){
		Point XC=centroid(X);
		DenseMatrix64F XM=new DenseMatrix64F(X.numRows, X.numCols);
		for(int i=0; i<X.numRows; i++){
			XM.set(i, 0, X.get(i, 0)-XC.x);
			XM.set(i, 1, X.get(i, 1)-XC.y);
			XM.set(i, 2, X.get(i, 2)-XC.z);
		}
		return XM;
	}
	
	public static DenseMatrix64F translation(DenseMatrix64F X, Point XC){
		DenseMatrix64F XM=new DenseMatrix64F(X.numRows, X.numCols);
		for(int i=0; i<X.numRows; i++){
			XM.set(i, 0, X.get(i, 0)-XC.x);
			XM.set(i, 1, X.get(i, 1)-XC.y);
			XM.set(i, 2, X.get(i, 2)-XC.z);
		}
		return XM;
	}
	
	public static double superimpose(DenseMatrix64F X, DenseMatrix64F Y){
		if(X.numRows != Y.numRows) throw new RuntimeException("From superimpose: the lengths of X and Y are not same");
		
		DenseMatrix64F XM=translation(X);
		DenseMatrix64F YM=translation(Y);

		DenseMatrix64F R=Kabsch(XM, YM);
		DenseMatrix64F XT=new DenseMatrix64F(XM.numCols, XM.numRows);
		CommonOps.transpose(XM, XT);
		DenseMatrix64F RX=new DenseMatrix64F(3, XT.numCols);
		CommonOps.mult(R, XT, RX);

		CommonOps.transpose(RX);
		return RMSD(RX, YM);
	}

	public static void print_R(DenseMatrix64F X, DenseMatrix64F Y) {
		if (X.numRows != Y.numRows) throw new RuntimeException("From superimpose: the lengths of X and Y are not same");

		DenseMatrix64F XM = translation(X);
		DenseMatrix64F YM = translation(Y);

		DenseMatrix64F R = Kabsch(XM, YM);
		System.out.print(R);
	}


	public static DenseMatrix64F List2Matrix(List<Point> coord){
		DenseMatrix64F ret=new DenseMatrix64F(coord.size(), 3);
		for(int i=0; i<coord.size(); i++){
			ret.set(i, 0, coord.get(i).x);
			ret.set(i, 1, coord.get(i).y);
			ret.set(i, 2, coord.get(i).z);
		}
		return ret;
	}

	public static double cal_rmsd_from_index(List<Integer> Map1_index, List<Integer> Map2_index, List<Point> coord1, List<Point> coord2){
		DenseMatrix64F Map1_X = new DenseMatrix64F(Map1_index.size(), 3);
		DenseMatrix64F Map2_Y = new DenseMatrix64F(Map2_index.size(), 3);

		for (int i = 0; i < Map1_index.size(); i++) {
			Map1_X.set(i, 0, coord1.get(Map1_index.get(i)).x);
			Map1_X.set(i, 1, coord1.get(Map1_index.get(i)).y);
			Map1_X.set(i, 2, coord1.get(Map1_index.get(i)).z);

			Map2_Y.set(i, 0, coord2.get(Map2_index.get(i)).x);
			Map2_Y.set(i, 1, coord2.get(Map2_index.get(i)).y);
			Map2_Y.set(i, 2, coord2.get(Map2_index.get(i)).z);
		}

		double rmsd = Geom.superimpose(Map1_X, Map2_Y);
		return rmsd;
	}

}
