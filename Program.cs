using System;
using System.Collections.Generic;
using System.Security.Cryptography;
using System.IO;
using System.Text;
using System.Drawing;

namespace SpatialHeterogeneityAnalysis
{
    class Program
    {
        static void Main(string[] args)
        {
            //Console.WriteLine("Hello World!");

            // CellArea から細胞間距離 l_c を計算
            double lc;
            if (GetValue("Average cell area (cm^2)", out double val))
            {
                lc = Math.Sqrt(2.0 / Math.Sqrt(3.0) * val) * 1E4;
                Console.WriteLine("Distance between cells l_c = " + lc + " (um)");
            }
            else return;
            Console.WriteLine();
            
            // 容器半径から Xsize, Ysize を計算
            int xSize, ySize;
            if (GetValue("Radius of culture vessel (mm)", out double radius))
            {
                double v = radius * 2 * 1E3 / lc;
                xSize = (int)(v + 2);
                ySize = (int)(v / Delta.Cf_y + 2);
                Console.WriteLine("Xsize = " + xSize + " (-)");
                Console.WriteLine("Ysize = " + ySize + " (-)");
            }
            else return;
            Console.WriteLine();

            // 播種密度を入力
            if (!GetValue("Seeding density (cells/cm^2)", out double density)) return;
            Console.WriteLine();

            // 画像保存用
            if (!GetValue("Image captur interval (mm)", out double interval)) return;
            Console.WriteLine();

            // 準備
            Common.Preparation();
            BoundaryConditions.SetParameter(false, false, false);
            CultureSpace.SetParameter(xSize, ySize, 3, lc, 1);
            // マップ作成
            if (!CultureSpace.MapCreation())
            { return; }
            // 細胞の配置
            if (!Seeding.Run_Random_CenterBiased(out List<CellData> cells, density))
            //if (!Seeding.Run(out List<CellData> cells))
            { return; }

            Capturing capturing = new Capturing();
            capturing.Run(interval);

            // 保存
            Output.Run(cells);
        }

        private static bool GetValue(string str, out double val, bool f = false)
        {
            Console.WriteLine("Set value: " + str);
            string s = Console.ReadLine();
            if (double.TryParse(s, out val))
            {
                return true;
            }
            else
            {
                Console.WriteLine("Error");
                if (f)
                {
                    return GetValue(str, out val, f);
                }
                else
                {
                    return false;
                }
            }
        }

        private static bool GetValue(string str, out int val, bool f = false)
        {
            Console.WriteLine("Set value: " + str);
            string s = Console.ReadLine();
            if (int.TryParse(s, out val))
            {
                return true;
            }
            else
            {
                Console.WriteLine("Error");
                if (f)
                {
                    return GetValue(str, out val, f);
                }
                else
                {
                    return false;
                }
            }
        }
    }
    public static class Common
    {
        public static void Preparation()
        {
            DirectionWeight();
            Rand = new Random();
        }

        // 静的自動プロパティ
        /// <summary>
        /// 方向の重み
        /// </summary>
        public static int[] Omega { get; private set; }
        private static Random Rand { get; set; } = new Random(); // 2019.12.19
        /// <summary>
        /// 0.0 以上 1.0 未満のランダムな浮動小数点数を返します。
        /// </summary>
        /// <returns></returns>
        public static double Rand_NextDouble() // 2019.12.19
        {
            lock (LockObj) { return Rand.NextDouble(); }
        }
        /// <summary>
        /// 指定した最大値より小さい 0 以上のランダムな整数を返します。
        /// </summary>
        /// <param name="maxValue"></param>
        /// <returns></returns>
        public static int Rand_Next(int maxValue) // 2020.02.26
        {
            lock (LockObj) { return Rand.Next(maxValue); }
        }
        private static readonly object LockObj = new(); // 2019.12.19

        private static void DirectionWeight()
        {
            //Omega = new int[20] {
            //1, 1, 1, 3, 1, 1, 1,
            //3, 3, 3,    3, 3, 3,
            //1, 1, 1, 3, 1, 1, 1 };

            //Omega = new int[20] { // 単位球内（内部の六角柱も含む）
            //     113,  113,  113, 884, 113,  113,  113,
            //    1146, 1146, 1146,     1146, 1146, 1146,
            //     113,  113,  113, 884, 113,  113,  113 };

            //Omega = new int[20] { // 単位球内（20近傍）
            //     143,  143,  143, 1112, 143,  143,  143,
            //    1010, 1010, 1010,      1010, 1010, 1010,
            //     143,  143,  143, 1112, 143,  143,  143 };

            Omega = new int[20] { // 単位球面
                292, 292, 292, 746, 292, 292, 292,
                834, 834, 834,      834, 834, 834,
                292, 292, 292, 746, 292, 292, 292 };
        }

        public static Delta RandomDirection(List<(Delta, int)> weight) // 2020.03.05
        {
            int[] arr = new int[weight.Count];
            int sum = 0;
            for (int i = 0; i < arr.Length; i++)
            {
                (_, int w) = weight[i];
                arr[i] = w;
                sum += arr[i];
            }
            if (sum > 0)
            {
                double val = sum * Rand_NextDouble();
                for (int i = 0; i < arr.Length; i++)
                {
                    (Delta d, int w) = weight[i];
                    val -= w;
                    if (val < 0)
                    { return d; }
                }
            }
            return null;
        }
        public static Delta RandomDirection(double[] weight)
        {
            double sum = 0;
            for (int i = 0; i < 20; i++)
            { sum += weight[i]; }
            if (sum > 0)
            {
                double val = sum * Rand_NextDouble();
                for (int i = 0; i < 20; i++)
                {
                    val -= weight[i];
                    if (val < 0)
                    { return Delta.GetDelta(i); }
                }
            }
            return null;
        }

        public static bool IsDeadCell(CellData c)
        {
            return c == null;
        }


        public static int[] RandomlySort(List<CellData> cells) // 2020.03.05
        {
            List<int> arr = new List<int>();
            for (int i = 0; i < cells.Count; i++)
            {
                if (!IsDeadCell(cells[i]))
                { arr.Add(cells[i].Index); }
            }
            for (int i = 0; i < arr.Count; i++)
            {
                int val = (int)(Rand_NextDouble() * (arr.Count - i)) + i;
                int tmp = arr[i];
                arr[i] = arr[val];
                arr[val] = tmp;
            }
            return arr.ToArray();
        }

        public static int[] RandomlySort(object[] obj)
        {
            int len = obj.Length;
            int[] arr = new int[len];
            for (int i = 0; i < len; i++)
            { arr[i] = i; }
            for (int i = 0; i < len; i++)
            {
                int val = (int)(Rand_NextDouble() * (len - i)) + i;
                int tmp = arr[i];
                arr[i] = arr[val];
                arr[val] = tmp;
            }
            return arr;
        }

        public static double UniformWhiteNoise(double rate)
        {
            // 範囲：[1-rate, 1+rate)
            return (Rand_NextDouble() * 2.0 - 1.0) * rate + 1.0;
        }

        public static int CellReplacement(List<CellData> cells, CellData c, Delta d)
        {
            Point3D cp = c.Location;
            Point3D p = c.Location + d;
            //Point3D p = BoundaryConditions.Check(c.Location, d);
            int ind2 = CultureSpace.GetMap(p);
            CultureSpace.SetMap(p, CultureSpace.GetMap(c.Location));
            CultureSpace.SetMap(c.Location, ind2);
            c.Location = p;
            if (ind2 >= 0) // 細胞なら
            {
                CellData c2 = CellData.Find(cells, ind2);
                c2.Location = cp;
            }
            return ind2;
        }
        ///// <summary>
        ///// ターゲット細胞が基質単位キューブとの細胞基質間結合エネルギーを持たないか
        ///// </summary>
        ///// <param name="c"></param>
        ///// <returns></returns>
        //public static bool IsDetatched_cs(CellData c) // 20191129
        //{
        //    int val = CultureSpace.GetMap(c.Location, new Delta(0, 0, -1));
        //    //int val = CultureSpace.MapGet(
        //    //    BoundaryConditions.Check(c.Location.X, c.Location.Y, c.Location.Z - 1));
        //    if (val == -3) // 注目細胞の直下が基質単位キューブならば
        //    {
        //        if (SubstrateAbility.Flag)
        //        {
        //            // 基質モジュールを使う場合
        //            for (int i = 13; i < 20; i++)
        //            {
        //                if (BasicConnectionEnergy.GetEcs(c, Delta.GetDelta(i)) > 0.0)
        //                {
        //                    // 基質との結合が1つでもあれば
        //                    return false;
        //                }
        //            }
        //            // 基質との結合がなければ
        //            return true;
        //        }
        //        else
        //        { return c.E_cs == 0.0; }
        //    }
        //    else
        //    {
        //        // 注目細胞の直下に基質単位キューブがなければ
        //        return true;
        //    }
        //}

        ///// <summary>
        ///// indexの細胞の近傍の細胞の細胞種を探索する。細胞種→0以上、基質→-1、それ以外→Int32最小値
        ///// </summary>
        ///// <param name="cells"></param>
        ///// <param name="index"></param>
        ///// <returns></returns>
        //public static int[] GetSurroundingCellType(List<CellData> cells, int index)
        //{
        //    CellData c = CellData.Find(cells, index);
        //    if (IsDeadCell(c)) { return null; }
        //    int[] result = new int[20];
        //    //for (int i = 0; i < 20; i++)
        //    _ = Parallel.For(0, 20, i =>
        //    {
        //        Point3D p = c.Location + Delta.GetDelta(i);
        //        //Point3D p = BoundaryConditions.Check(c.Location, Delta.GetDelta(i));
        //        int val = CultureSpace.GetMap(p);
        //        if (val >= 0)
        //        { result[i] = CellData.Find(cells, val).Cell_T; }
        //        else if (val == -3)
        //        { result[i] = -1; }
        //        else
        //        { result[i] = int.MinValue; }
        //    });
        //    return result;
        //}

        //public static void InitializeParameter<T>(int typeNum, List<bool> flag, T defVal, List<T> param)// where T : IComparable
        //{
        //    if (param == null)
        //    {
        //        param = new List<T>();
        //        flag = new List<bool>();
        //    }
        //    while (param.Count <= typeNum)
        //    {
        //        param.Add(defVal);
        //        if (param.Count > flag.Count)
        //        { flag.Add(false); }
        //    }
        //}
        //public static void InitializeParameter<T>(int typeNum, List<bool> flag, List<List<T>> param)// where T : IComparable
        //{
        //    if (param == null)
        //    {
        //        param = new List<List<T>>();
        //        flag = new List<bool>();
        //    }
        //    while (param.Count <= typeNum)
        //    {
        //        param.Add(new List<T>());
        //        if (param.Count > flag.Count)
        //        { flag.Add(false); }
        //    }
        //}
        //public static void InitializeParameter<T>(int typeNum, T defVal, params List<T>[] param)// where T : IComparable
        //{
        //    for (int i = 0; i < param.Length; i++)
        //    {
        //        if (param[i] == null)
        //        {
        //            param[i] = new List<T>();
        //        }
        //        while (param[i].Count <= typeNum)
        //        {
        //            param[i].Add(defVal);
        //        }
        //    }
        //}
    }
    public class CellData
    {
        /// <summary>
        /// CellData の新規作成
        /// </summary>
        public CellData(sbyte cellT = 0) // 2020.12.03
        {
            // 細胞インデックス
            Index = -1;
            // 細胞位置
            Location = new Point3D(-1, -1, -1);
        }

        /// <summary>
        /// CellData c のディープコピー
        /// </summary>
        /// <param name="c">コピー元のCellData</param>
        public CellData(CellData c)
        {
            // 細胞インデックス
            Index = c.Index;
            // 細胞位置
            Location = new Point3D(c.Location);//.X, c.Location.Y, c.Location.Z);
        }

        /// <summary>
        /// 細胞のインデックスを取得または設定します。
        /// </summary>
        internal int Index { get; set; }
        /// <summary>
        /// 細胞の座標を取得または設定します。
        /// </summary>
        internal Point3D Location { get; set; }

        internal CellData Clone()
        {
            return Clone(this);
        }

        internal static CellData Clone(CellData c)
        {
            return c == null ? null : new CellData(c);
        }

        /// <summary>
        /// cellsのindexに対応する配列のインデックスを出力する
        /// </summary>
        /// <param name="cells"></param>
        /// <param name="index"></param>
        /// <returns></returns>
        internal static CellData Find(List<CellData> cells, int index)
        {
            return cells[index];
        }

        internal static int FindIndex(List<CellData> cells, int index)
        {
            return index;
        }

    }
    public class Point3D
    {
        /// <summary>
        /// 新しいインスタンスを初期化、 CellData.Point 座標を指定しています。
        /// </summary>
        public Point3D(Point3D p) { X = p.X; Y = p.Y; Z = p.Z; }
        /// <summary>
        /// 新しいインスタンスを初期化、 CellData.Point 座標を指定しています。
        /// </summary>
        /// <param name="x">ポイントの横方向の位置。</param>
        /// <param name="y">ポイントの縦方向の位置。</param>
        /// <param name="z">ポイントの高さ方向の位置。</param>
        public Point3D(int x, int y, int z) { X = x; Y = y; Z = z; }
        /// <summary>
        /// Point3D の x 座標を取得または設定します。
        /// </summary>
        public int X { get; set; }
        /// <summary>
        /// Point3D の y 座標を取得または設定します。
        /// </summary>
        public int Y { get; set; }
        /// <summary>
        /// Point3D の z 座標を取得または設定します。
        /// </summary>
        public int Z { get; set; }

        /// <summary>
        /// DeltaからCellData.Pointにキャストします。
        /// </summary>
        /// <param name="v"></param>
        public static explicit operator Point3D(Delta v)
        {
            return new Point3D(v.DX, v.DY, v.DZ);
        }
        public static bool IsSameLocation(Point3D p1, Point3D p2)
        {
            if (p1.X == p2.X && p1.Y == p2.Y && p1.Z == p2.Z)
            { return true; }
            else
            { return false; }
        }
        public bool IsSameLocation(Point3D p)
        {
            if (X == p.X && Y == p.Y && Z == p.Z)
            { return true; }
            else
            { return false; }
        }
        public static Point3D GetPointWithCheckBoundaryConditions(int x, int y, int z)
        {
            return BoundaryConditions.Check(x, y, z);
        }
        public static Point3D GetPointWithCheckBoundaryConditions(Point3D point, Delta delta)
        {
            return point + delta;
        }
        public static Point3D operator +(Point3D left, Delta right)
        {
            return BoundaryConditions.Check(left, right);
        }
        public static Point3D operator +(Delta left, Point3D right)
        {
            return BoundaryConditions.Check(right, left);
        }

    }
    public class Delta
    {
        public Delta(int dx, int dy, int dz)
        { DX = dx; DY = dy; DZ = dz; }
        public Delta()
        { DX = 0; DY = 0; DZ = 0; }

        public int DX { get; set; }
        public int DY { get; set; }
        public int DZ { get; set; }
        /// <summary>
        /// 六角格子座標から実座標へのy軸の補正係数Correction factor: sqrt(3)/2
        /// </summary>
        public static double Cf_y { get; } = Math.Sqrt(3) / 2.0;

        public static Delta Copy(Delta src)
        {
            return new Delta(src.DX, src.DY, src.DZ);
        }
        public static void Copy(Delta src, Delta dsc)
        {
            dsc.DX = src.DX;
            dsc.DY = src.DY;
            dsc.DZ = src.DZ;
        }

        public static Delta GetDelta(Direction.DIR dir)
        {
            switch (dir)
            {
                case Direction.DIR.UL2: return new Delta(-1, -1, 1);
                case Direction.DIR.UR2: return new Delta(1, -1, 1);
                case Direction.DIR.L_2: return new Delta(-2, 0, 1);
                case Direction.DIR.C_2: return new Delta(0, 0, 1);
                case Direction.DIR.R_2: return new Delta(2, 0, 1);
                case Direction.DIR.LL2: return new Delta(-1, 1, 1);
                case Direction.DIR.LR2: return new Delta(1, 1, 1);
                case Direction.DIR.UL1: return new Delta(-1, -1, 0);
                case Direction.DIR.UR1: return new Delta(1, -1, 0);
                case Direction.DIR.L_1: return new Delta(-2, 0, 0);
                case Direction.DIR.C_1: return new Delta(0, 0, 0);
                case Direction.DIR.R_1: return new Delta(2, 0, 0);
                case Direction.DIR.LL1: return new Delta(-1, 1, 0);
                case Direction.DIR.LR1: return new Delta(1, 1, 0);
                case Direction.DIR.UL0: return new Delta(-1, -1, -1);
                case Direction.DIR.UR0: return new Delta(1, -1, -1);
                case Direction.DIR.L_0: return new Delta(-2, 0, -1);
                case Direction.DIR.C_0: return new Delta(0, 0, -1);
                case Direction.DIR.R_0: return new Delta(2, 0, -1);
                case Direction.DIR.LL0: return new Delta(-1, 1, -1);
                case Direction.DIR.LR0: return new Delta(1, 1, -1);
                default: return null;
            }
        }
        public static Delta GetDelta(int i)
        {
            return GetDelta((Direction.DIR)i);
        }
        public static Delta GetDelta(Point3D from, Point3D to) // 2021.08.04
        {
            return new Delta(to.X - from.X, to.Y - from.Y, to.Z - from.Z);
        }

        /// <summary>
        /// dだけ離れた2点間の距離の2乗を計算する。単位は無次元（µm^2でない）。
        /// </summary>
        /// <returns></returns>
        public static double GetLength_pow2(Delta d)
        {
            return d.DX * d.DX / 4.0 + d.DY * d.DY * Cf_y * Cf_y + d.DZ * d.DZ;
        }
        /// <summary>
        /// dだけ離れた2点間のグリッド距離を計算する。単位は無次元（µmでない）。
        /// </summary>
        /// <returns></returns>
        public static double GetLength(Delta d)
        {
            return Math.Sqrt(GetLength_pow2(d));
        }
        /// <summary>
        /// 原点からの線分と点との最短距離の2乗
        /// </summary>
        /// <param name="distance"></param>
        /// <param name="point"></param>
        /// <returns></returns>
        public static double DistanceFromPointToLine(Delta distance, Delta point)
        {
            double xx = distance.DX * distance.DX / 4.0;
            double yy = distance.DY * distance.DY * Cf_y * Cf_y;
#if false
            double r2 = xx + yy;
            double len2; // 任意の点と任意の線分との最短距離の2乗
            double tt = distance.DX * point.DX / 4.0 + distance.DY * point.DY * Cf_y * Cf_y;
            if (tt < 0)
            { len2 = GetLength_pow2(point); }
            else if (tt > r2)
            { len2 = GetLength_pow2(distance - point); }
            else
            {
                double f1 = ((distance.DY * point.DX) - (distance.DX * point.DY)) / 2.0 * Cf_y;
                len2 = f1 * f1 / r2;
            }
            return len2;
#elif false
            double r2 = xx + yy;
            double f1 = ((distance.DY * point.DX) - (distance.DX * point.DY)) / 2.0 * Cf_y;
            return r2 == 0 ? 0 : f1 * f1 / r2;
#else
            double zz = distance.DZ * distance.DZ;
            double dxpy = distance.DX * point.DY / 2.0 * Cf_y;
            double dxpz = distance.DX * point.DZ / 2.0;
            double dypx = distance.DY * point.DX / 2.0 * Cf_y;
            double dypz = distance.DY * point.DZ * Cf_y;
            double dzpx = distance.DZ * point.DX / 2.0;
            double dzpy = distance.DZ * point.DY * Cf_y;
            double xy = dxpy - dypx;
            double yz = dypz - dzpy;
            double zx = dzpx - dxpz;
            double d = xy * xy + yz * yz + zx * zx;
            return d / (xx + yy + zz);
#endif
        }

        #region operator
        public static Delta operator -(Delta left, Delta right)
        {
            return new Delta(left.DX - right.DX, left.DY - right.DY, left.DZ - right.DZ);
        }
        public static Delta operator +(Delta left, Delta right)
        {
            return new Delta(left.DX + right.DX, left.DY + right.DY, left.DZ + right.DZ);
        }
        public static Delta operator -(Delta d)
        {
            return new Delta(-d.DX, -d.DY, -d.DZ);
        }
        /// <summary>
        /// Point3D から Delta にキャストします。
        /// </summary>
        /// <param name="p">Location</param>
        public static explicit operator Delta(Point3D p)
        {
            return new Delta(p.X, p.Y, p.Z);
        }
        #endregion

    }
    public class Direction
    {
        /// <summary>
        /// upper layer    middle layer    under layer
        ///   0   1          7   8           13  14
        /// 2   3   4      9  -1  10       15  16  17
        ///   5   6         11  12           18  19
        /// </summary>
        public enum DIR
        {
            /// <summary> upper left (upper layer) </summary>
            UL2 = 0,
            /// <summary> upper right (upper layer) </summary>
            UR2 = 1,
            /// <summary> left (upper layer) </summary>
            L_2 = 2,
            /// <summary> center (upper layer) </summary>
            C_2 = 3,
            /// <summary> right (upper layer) </summary>
            R_2 = 4,
            /// <summary> lower left (upper layer) </summary>
            LL2 = 5,
            /// <summary> lower right (upper layer) </summary>
            LR2 = 6,
            /// <summary> upper left (middle layer) </summary>
            UL1 = 7,
            /// <summary> upper right (middle layer) </summary>
            UR1 = 8,
            /// <summary> left (middle layer) </summary>
            L_1 = 9,
            /// <summary> center (middle layer) </summary>
            C_1 = -1,
            /// <summary> right (middle layer) </summary>
            R_1 = 10,
            /// <summary> lower left (middle layer) </summary>
            LL1 = 11,
            /// <summary> lower right (middle layer) </summary>
            LR1 = 12,
            /// <summary> upper left (underlayer) </summary>
            UL0 = 13,
            /// <summary> upper right (underlayer) </summary>
            UR0 = 14,
            /// <summary> left (underlayer) </summary>
            L_0 = 15,
            /// <summary> center (under layer) </summary>
            C_0 = 16,
            /// <summary> right (underlayer) </summary>
            R_0 = 17,
            /// <summary> lower left (underlayer) </summary>
            LL0 = 18,
            /// <summary> lower right (underlayer) </summary>
            LR0 = 19,
            /// <summary> No direction </summary>
            NULL = int.MinValue,
        }

        public int[] DirectionToVariation(DIR d)
        { // [x, y, z]
            int[] v;
            switch (d)
            {
                case DIR.UL2: v = new int[] { -1, -1, 1 }; break;
                case DIR.UR2: v = new int[] { 1, -1, 1 }; break;
                case DIR.L_2: v = new int[] { -2, 0, 1 }; break;
                case DIR.C_2: v = new int[] { 0, 0, 1 }; break;
                case DIR.R_2: v = new int[] { 2, 0, 1 }; break;
                case DIR.LL2: v = new int[] { -1, 1, 1 }; break;
                case DIR.LR2: v = new int[] { 1, 1, 1 }; break;
                case DIR.UL1: v = new int[] { -1, -1, 0 }; break;
                case DIR.UR1: v = new int[] { 1, -1, 0 }; break;
                case DIR.L_1: v = new int[] { -2, 0, 0 }; break;
                case DIR.R_1: v = new int[] { 2, 0, 0 }; break;
                case DIR.LL1: v = new int[] { -1, 1, 0 }; break;
                case DIR.LR1: v = new int[] { 1, 1, 0 }; break;
                case DIR.UL0: v = new int[] { -1, -1, -1 }; break;
                case DIR.UR0: v = new int[] { 1, -1, -1 }; break;
                case DIR.L_0: v = new int[] { -2, 0, -1 }; break;
                case DIR.C_0: v = new int[] { 0, 0, -1 }; break;
                case DIR.R_0: v = new int[] { 2, 0, -1 }; break;
                case DIR.LL0: v = new int[] { -1, 1, -1 }; break;
                case DIR.LR0: v = new int[] { 1, 1, -1 }; break;
                default: v = null; break;
            }
            return v;
        }
        public static DIR VariationToDirection(int[] arr)
        { // [x, y, z]
            if (arr[0] == -2)
            {
                if (arr[2] == -1) { return DIR.L_0; }
                else if (arr[2] == 0) { return DIR.L_1; }
                else if (arr[2] == 1) { return DIR.L_2; }
            }
            else if (arr[0] == -1)
            {
                if (arr[1] == -1)
                {
                    if (arr[2] == -1) { return DIR.UL0; }
                    else if (arr[2] == 0) { return DIR.UL1; }
                    else if (arr[2] == 1) { return DIR.UL2; }
                }
                else if (arr[1] == 1)
                {
                    if (arr[2] == -1) { return DIR.LL0; }
                    else if (arr[2] == 0) { return DIR.LL1; }
                    else if (arr[2] == 1) { return DIR.LL2; }
                }
            }
            else if (arr[0] == 0)
            {
                if (arr[1] == 0)
                {
                    if (arr[2] == -1) { return DIR.C_0; }
                    else if (arr[2] == 0) { return DIR.C_1; }
                    else if (arr[2] == 1) { return DIR.C_2; }
                }
                else
                {
                    return DIR.NULL;
                }
            }
            else if (arr[0] == 1)
            {
                if (arr[1] == -1)
                {
                    if (arr[2] == -1) { return DIR.UR0; }
                    else if (arr[2] == 0) { return DIR.UR1; }
                    else if (arr[2] == 1) { return DIR.UR2; }
                }
                else if (arr[1] == 1)
                {
                    if (arr[2] == -1) { return DIR.LR0; }
                    else if (arr[2] == 0) { return DIR.LR1; }
                    else if (arr[2] == 1) { return DIR.LR2; }
                }
            }
            else if (arr[0] == 2)
            {
                if (arr[2] == -1) { return DIR.R_0; }
                else if (arr[2] == 0) { return DIR.R_1; }
                else if (arr[2] == 1) { return DIR.R_2; }
            }
            return DIR.NULL;
        }

        public static DIR RandomDirection(double[,,] dr, Random rand)
        {
            double val = 0;
            for (int k = 0; k < dr.GetLength(0); k++)
            {
                for (int j = 0; j < dr.GetLength(1); j++)
                {
                    for (int i = 0; i < dr.GetLength(2); i++)
                    {
                        val += dr[k, j, i];
                    }
                }
            }
            // 0 <= val' < val
            while ((val *= rand.NextDouble()) == 1) ;

            for (int k = 0; k < dr.GetLength(0); k++)
            {
                for (int j = 0; j < dr.GetLength(1); j++)
                {
                    for (int i = 0; i < dr.GetLength(2); i++)
                    {
                        if ((val -= dr[k, j, i]) < 0)
                        { return VariationToDirection(new int[] { i - 1, j - 1, k - 1 }); }
                    }
                }
            }
            return DIR.NULL;
        }
        public static DIR RandomDirection(double[] dr, Random rand)
        {
            double val = 0;
            for (int i = 0; i < dr.Length; i++)
            {
                val += dr[i];
            }
            // 0 <= val' < val
            while ((val *= rand.NextDouble()) == 1) ;

            for (int i = 0; i < dr.Length; i++)
            {
                if ((val -= dr[i]) < 0)
                { return (DIR)i; }
            }
            return DIR.NULL;
        }

        public static DIR RandomDirection(Random rand)
        {
            double[] dr = new double[20];
            for (int i = 0; i < dr.Length; i++)
            {
                dr[i] = Common.Omega[i];
            }
            return RandomDirection(dr, rand);
        }

        public static DIR GetDirection(Delta d)
        {
            if (d.DZ == 1)
            {
                if (d.DY == -1)
                {
                    if (d.DX == -1) { return DIR.UL2; }
                    else if (d.DX == 1) { return DIR.UR2; }
                    else { return DIR.NULL; }
                }
                else if (d.DY == 0)
                {
                    if (d.DX == -2) { return DIR.L_2; }
                    else if (d.DX == 0) { return DIR.C_2; }
                    else if (d.DX == 2) { return DIR.R_2; }
                    else { return DIR.NULL; }
                }
                else if (d.DY == 1)
                {
                    if (d.DX == -1) { return DIR.LL2; }
                    else if (d.DX == 1) { return DIR.LR2; }
                    else { return DIR.NULL; }
                }
                else { return DIR.NULL; }
            }
            else if (d.DZ == 0)
            {
                if (d.DY == -1)
                {
                    if (d.DX == -1) { return DIR.UL1; }
                    else if (d.DX == 1) { return DIR.UR1; }
                    else { return DIR.NULL; }
                }
                else if (d.DY == 0)
                {
                    if (d.DX == -2) { return DIR.L_1; }
                    else if (d.DX == 0) { return DIR.C_1; }
                    else if (d.DX == 2) { return DIR.R_1; }
                    else { return DIR.NULL; }
                }
                else if (d.DY == 1)
                {
                    if (d.DX == -1) { return DIR.LL1; }
                    else if (d.DX == 1) { return DIR.LR1; }
                    else { return DIR.NULL; }
                }
                else { return DIR.NULL; }
            }
            else if (d.DZ == -1)
            {
                if (d.DY == -1)
                {
                    if (d.DX == -1) { return DIR.UL0; }
                    else if (d.DX == 1) { return DIR.UR0; }
                    else { return DIR.NULL; }
                }
                else if (d.DY == 0)
                {
                    if (d.DX == -2) { return DIR.L_0; }
                    else if (d.DX == 0) { return DIR.C_0; }
                    else if (d.DX == 2) { return DIR.R_0; }
                    else { return DIR.NULL; }
                }
                else if (d.DY == 1)
                {
                    if (d.DX == -1) { return DIR.LL0; }
                    else if (d.DX == 1) { return DIR.LR0; }
                    else { return DIR.NULL; }
                }
                else { return DIR.NULL; }
            }
            else { return DIR.NULL; }
        }

        public static DIR ReverseDirection(DIR d)
        {
            return d == (DIR)(-1) ? (DIR)(-1) : 19 - d;
        }
    }
    public static class BoundaryConditions
    {
        public static void SetParameter(bool periodic_X, bool periodic_Y, bool periodic_Z)
        {
            PB_X = periodic_X;
            PB_Y = periodic_Y;
            PB_Z = periodic_Z;
        }

        // 静的プロパティ
        private static bool PB_X { get; set; }
        private static bool PB_Y { get; set; }
        private static bool PB_Z { get; set; }

        /// <summary>
        /// 周期境界条件による座標の補正
        /// </summary>
        /// <param name="val">補正前の座標</param>
        /// <param name="pb">周期境界条件</param>
        /// <param name="arrSize">用意する計算空間</param>
        /// <returns>補正後の座標，+/-方向に境界を超えた回数</returns>
        private static (int, int) Check(int val, bool pb, int arrSize) // 2020.11.04
        {
            if (pb)
            {
                arrSize--;
                arrSize--;
                //int res = val % arrSize;
                int p = (int)Math.Floor((double)val / arrSize); // ゼロでなければ周期境界になっている(periodicity)
                int res = val - arrSize * p;
                if (res == 0) { return (arrSize, p - 1); }
                else { return (res, p); }
            }
            else
            {
                arrSize--;
                if (val > 0)
                {
                    if (arrSize > val) { return (val, 0); }
                    else { return (arrSize, 0); }
                }
                else { return (0, 0); }
            }
        }

        /// <summary>
        /// Boundary Condition Checker
        /// </summary>
        /// <param name="x">X座標</param>
        /// <param name="y">Y座標</param>
        /// <param name="z">Z座標</param>
        /// <returns></returns>
        internal static Point3D Check(int x, int y, int z)
        {
            (int outy, int py) = Check(y, PB_Y, CultureSpace.Ysize);
            int x0 = 0; // y軸が周期境界のときのx軸の補正係数
            // y軸が周期境界かつYsizeが奇数なら // 2020.11.04
            if (py != 0 && CultureSpace.Ysize % 2 == 1)
            {
                // 壁面キューブがあるからyの取りうる範囲は[1, Ysize - 2]。それ以外ならxを補正する。
                if (y <= 0 || y > CultureSpace.Ysize - 1)
                { x0 = py; }
            }
            (int outx, _) = Check((int)Math.Floor((x + x0) / 2.0), PB_X, CultureSpace.Xsize); // 2020.10.30
            (int outz, _) = Check(z, PB_Z, CultureSpace.Zsize);
            // yの補正値が奇数か偶数かでxの値を補正する
            if (outy % 2 == 0)
            { return new Point3D(outx * 2, outy, outz); }
            else
            { return new Point3D(outx * 2 + 1, outy, outz); }
        }

        /// <summary>
        /// Boundary Condition Checker
        /// </summary>
        /// <param name="point">座標</param>
        /// <param name="delta">移動量</param>
        /// <returns></returns>
        internal static Point3D Check(Point3D point, Delta delta)
        {
            int x = point.X + delta.DX;
            int y = point.Y + delta.DY;
            int z = point.Z + delta.DZ;
            return Check(x, y, z);
        }
    }
    public static class CultureSpace
    {
        public static void SetParameter(int xSize, int ySize, int zSize, double lc, int type = 1)
        {
            Size_lc = lc;
            Xsize = xSize;
            Ysize = ySize;
            Zsize = zSize;
            _type = (TYPE)type;
        }

        // 静的メンバー変数
        //private static string _mapType;
        private static TYPE _type;
        private static int[,,] Map;

        // 静的プロパティ
        internal static double Size_lc { get; private set; }
        //internal static double Size_hc { get; private set; }
        internal static int Xsize { get; private set; }
        internal static int Ysize { get; private set; }
        internal static int Zsize { get; private set; }

        /// <summary>
        /// マップが有効かどうかを示します。
        /// </summary>
        internal static bool MapEnabled => Map != null;

        private enum TYPE
        {
            Square = 0,
            Circle = 1,
        }

        public static bool MapCreation()
        {
            switch (_type)
            {
                case TYPE.Square:
                    Creation_Square();
                    return true;
                case TYPE.Circle:
                    Creation_Circle();
                    return true;
                default:
                    return false;
            }
        }

        private static void Creation_Square()
        {
            Map = new int[Zsize, Ysize, Xsize];

            for (int j = 1; j < Ysize - 1; j++)
            {
                for (int i = 1; i < Xsize - 1; i++)
                {
                    SetMap(i * 2, j, Zsize - 1, -2); // block
                }
            }
            for (int k = 1; k < Zsize - 1; k++)
            {
                for (int j = 1; j < Ysize - 1; j++)
                {
                    for (int i = 1; i < Xsize - 1; i++)
                    {
                        SetMap(i * 2, j, k, -1); // empty grid
                    }
                }
            }
            for (int k = 0; k < Zsize; k++)
            {
                for (int j = 0; j < Ysize; j++)
                {
                    SetMap(0, j, k, -2); // block
                    SetMap((Xsize - 1) * 2, j, k, -2); // block
                }
            }
            for (int k = 0; k < Zsize; k++)
            {
                for (int i = 0; i < Xsize; i++)
                {
                    SetMap(i * 2, 0, k, -2); // block
                    SetMap(i * 2, Ysize - 1, k, -2); // block
                }
            }
            for (int j = 0; j < Ysize; j++)
            {
                for (int i = 0; i < Xsize; i++)
                {
                    SetMap(i * 2, j, 0, -3); // substrate
                }
            }
        }
        private static void Creation_Circle()
        {
            Map = new int[Zsize, Ysize, Xsize];

            int y2 = Ysize / 2 - 1;
            double x2 = y2 % 2 == 0 ? Xsize / 2 - 1 : Xsize / 2 - 0.5;
            double y3 = Ysize * Delta.Cf_y / 2.0;
            double r = (y3 < x2 ? y3 : x2) - 1; // 半径はXとYで小さい方に合わせる

            for (int j = 0; j < Ysize; j++)
            {
                for (int i = 0; i < Xsize; i++)
                {
                    double len = j % 2 == 0
                        ? Delta.GetLength(new Delta((int)((i - x2) * 2), j - y2, 0))
                        : Delta.GetLength(new Delta((int)((i - x2) * 2 + 1), j - y2, 0));
                    if (len <= r)
                    {
                        SetMap(i * 2, j, 0, -3);
                        for (int k = 1; k < Zsize - 1; k++)
                        {
                            SetMap(i * 2, j, k, -1);
                        }
                    }
                    else
                    {
                        SetMap(i * 2, j, 0, -2);
                        SetMap(i * 2, j, 1, -3);
                    }
                    SetMap(i * 2, j, Zsize - 1, -2);
                }
            }
        }

        internal static int GetMap(int row, int col, int dep)
        {
            return Map[dep, col, row / 2];
        }
        internal static int GetMap(Point3D point)
        {
            return GetMap(point.X, point.Y, point.Z);
        }
        internal static int GetMap(Point3D point, Delta delta) // 2020.03.04
        {
            //return MapGet(point.X + delta.DX, point.Y + delta.DY, point.Z + delta.DZ);
            return GetMap(point + delta);
        }

        internal static void SetMap(int row, int col, int dep, int value)
        {
            Map[dep, col, row / 2] = value;
        }
        internal static void SetMap(Point3D point, int value)
        {
            SetMap(point.X, point.Y, point.Z, value);
        }

        internal static bool IsCorrect(Point3D point, int ind, out int result)
        {
            result = GetMap(point);
            return result == ind;
        }
        internal static bool IsCorrect(List<CellData> cells, int id, out int result)
        {
            return IsCorrect(cells[id].Location, cells[id].Index, out result);
        }
    }
    public class Seeding
    {
        //private enum InitializationMethod
        //{
        //    CSV = 0, // Input CellsData
        //    NEW = 1, // Initialize CellsData
        //}
        private enum SeedingTYPE
        {
            Random_NoBias = 0,
            Random_CenterBiased = 1,
            Uniform = 2,
            SingleColony_Layer = 3,
            SingleColony_CellNumber = 4,
            MultipleColony_Random_mean_colony_size = 5, // 2019.12.25
            MultupleColony_Random_inoculum_size = 6, // 2019.12.25
            MultipleColony_Determined = 7,
            MultipleAggregation_Random = 8,
            LoadCellPlacement = 9,
            NonUniformDistribution = 10, // 2021.07.07
        }

        // メンバー変数
        //private static InitializationMethod initMet;
        //private static string InputCSV;
        //private static string InputCSVpre; // 2020.05.26
        private static SeedingTYPE sType = SeedingTYPE.Random_NoBias;
        private static List<double> X_0; // Inoculum size 接種細胞密度
        private static int ColonyRadius = 1; // コリニーの半径
        private static int ColonyCells = 1; // コロニー内の細胞数
        private static int GridInterval = 1; // Grid interval
        private static int NumOfColonies = 1;
        private static int ColonyNum_SD = 0; // 2020.03.11
        private static int ColonyCellsMu = 1;
        private static int ColonyCellsSigma = 0;
        //private static string LoadGeomInfo;
        //private static string ImagePath;
        //private static double PixelSize; // 2020.05.27
        private static double Uniformity; // 2021.07.07
        private static bool Uniformity_Inf; // 2021.07.07
        //private static bool ErrFlag;
        private static double ConcentrationParameter; // 2021.08.02
        private static double CenterPeripheryRatio; // 2021.08.31

        // 実処理
        public static bool Run(out List<CellData> cells)
        {
            //bool preData = false; // 2020.05.26
            //switch (initMet)
            //{
            //    case InitializationMethod.CSV:
            //        (cells, preData) = InputCellsData(); // 2020.05.26
            //        break;
            //    case InitializationMethod.NEW:
            switch (sType)
            {
                case SeedingTYPE.Random_NoBias:
                    cells = CellInitialization_Random();
                    if (cells == null) { break; }
                    InitialCellPlacement_Random(cells);
                    break;
                case SeedingTYPE.Random_CenterBiased:
                    cells = CellInitialization_Random();
                    if (cells == null) { break; }
                    InitialCellPlacement_Random_CenterBiased(cells);
                    break;
                case SeedingTYPE.Uniform:
                    cells = CellInitialization_Uniform();
                    break;
                case SeedingTYPE.SingleColony_Layer:
                    cells = CellInitialization_SingleColony();
                    InitialCellPlacement_SingleColony(cells);
                    break;
                case SeedingTYPE.SingleColony_CellNumber:
                    cells = CellInitialization_SingleColony_CellNumber();
                    InitialCellPlacement_SingleColony_CellNumber(cells);
                    break;
                case SeedingTYPE.MultipleColony_Random_mean_colony_size:
                    cells = CellInitialization_MultipleColony_Random_MeanColony();
                    break;
                case SeedingTYPE.MultupleColony_Random_inoculum_size: // 2019.12.25
                    cells = CellInitialization_MultipleColony_Random_Inoculum();
                    break;
                case SeedingTYPE.MultipleAggregation_Random:
                    cells = CellInitialization_MultipleAggregation_Random();
                    break;
                //case SeedingTYPE.MultipleColony_Determined:
                //    Console.WriteLine("未実装：MultipleColony_Determined");
                //    cells = null;
                //    break;
                //case SeedingTYPE.LoadCellPlacement:
                //    cells = CellInitialization_LoadCellPlacement(); // 2020.05.27
                //    break;
                case SeedingTYPE.NonUniformDistribution: // 2021.07.07
                    cells = CellInitialization_Random();
                    if (cells == null) { break; }
                    InitialCellPlacement_NonUniform(ref cells);
                    break;
                default:
                    cells = null;
                    break;
            }
            //        break;
            //    default:
            //        cells = null;
            //        break;
            //}

            //if (cells != null && !preData)
            //{
            //    // 初期状態のCellDataとMapは現在の状態のデータを使用する
            //    Common.MaintainStates(cells);
            //}

            return cells != null;
        }
        public static bool Run_Random_CenterBiased(out List<CellData> cells, double X0, double concentrationParameter = 0.5, double centerPeripheryRatio = 0.0)
        {
            X_0 = new List<double> { X0 };
            ConcentrationParameter = concentrationParameter;
            CenterPeripheryRatio = centerPeripheryRatio;

            cells = CellInitialization_Random();
            if (cells == null) { return false; }
            InitialCellPlacement_Random(cells);
            return cells != null;
        }
        public static bool Run_NonUniformDistribution(out List<CellData> cells, int uniformity, bool uniformity_Inf = false)
        {
            Uniformity = uniformity;
            Uniformity_Inf = uniformity_Inf;

            cells = CellInitialization_Random();
            if (cells == null) { return false; }
            InitialCellPlacement_NonUniform(ref cells);
            return cells != null;
        }

        //#region InputCellsData
        //private static (List<CellData>, bool) InputCellsData() // 2020.05.26
        //{
        //    bool flag = false;
        //    List<CellData> pre = InputCellsData(InputCSVpre);
        //    if (pre != null)
        //    {
        //        // 現在の状態のCellDataとMapを保持する
        //        Common.MaintainStates(pre);
        //        flag = true;
        //        // マップの初期化 2020.11.26
        //        CultureSpace.MapCreation();
        //    }

        //    return (InputCellsData(InputCSV), flag);
        //}
        //public static List<CellData> InputCellsData(string filename)
        //{
        //    if (!File.Exists(filename))
        //    { return null; }
        //    using (FileStream fs = new FileStream(filename, FileMode.Open, FileAccess.Read, FileShare.ReadWrite))
        //    using (StreamReader sr = new StreamReader(fs))
        //    {
        //        string[] str = new string[] {
        //            Output.Name_Index,
        //            Output.Name_X,
        //            Output.Name_Y,
        //            Output.Name_Z,
        //            Output.Name_hexX,
        //            Output.Name_hexY,
        //            Output.Name_hexZ,
        //            Output.Name_Cell_N,
        //            Output.Name_Cell_S,
        //            Output.Name_Cell_T,
        //            Output.Name_t_age,
        //            Output.Name_t_d,
        //            Output.Name_t_m,
        //            Output.Name_V_m,
        //            Output.Name_E_AJ,
        //            Output.Name_E_TJ,
        //            Output.Name_E_cs,
        //            Output.Name_E_max, // 2021.02.18
        //            Output.Name_Dir_am,
        //            Output.Name_Dir_pm,
        //            Output.Name_Dir_pd,
        //            Output.Name_THETA,
        //            Output.Name_T_de_act,
        //            Output.Name_T_de_total,
        //            Output.Name_CellMovement_len_app, // 2020.03.09
        //            Output.Name_CellMovement_len_act, // 2020.03.09
        //            Output.Name_CellMovement_len_pas, // 2020.03.09
        //            Output.Name_CellMovement_len_div, // 2020.03.09
        //            Output.Name_CellMovement_t_act, // 2020.03.09
        //            Output.Name_CellMovement_t_pas, // 2020.03.09
        //            Output.Name_CellMovement_t_div, // 2020.03.09
        //        };
        //        string[] labels = sr.ReadLine().Split(',');
        //        int[] lut = new int[str.Length];
        //        string errStr = "";
        //        for (int i = 0; i < str.Length; i++)
        //        {
        //            lut[i] = Array.FindIndex(labels, s => s == str[i]);
        //            // 入力されたCellsDataの列名が定義と異なる場合 // 20191204
        //            if (lut[i] == -1)
        //            { errStr += str[i] + "\r\n"; }
        //        }
        //        // 一度だけ確認する
        //        if (errStr != "" && ErrFlag)
        //        {
        //            if (MessageBox.Show(
        //                "Ignore the data of CellsData with the following column name?\r\n" + errStr,
        //                "Confirmation of cell seeding", MessageBoxButtons.YesNo, MessageBoxIcon.Exclamation)
        //                == DialogResult.Yes)
        //            { ErrFlag = false; }
        //            else
        //            { return null; }
        //        }
        //        if (labels[0] == Output.Name_Index)
        //        {
        //            List<CellData> cells = new List<CellData>();
        //            while (!sr.EndOfStream)
        //            {
        //                string[] strs = sr.ReadLine().Split(',');
        //                // CellTypeはここで定義 2020.12.03
        //                sbyte cellT = (sbyte)(lut[9] > -1 ? Input_SByte(strs[lut[9]]) : 0);
        //                CellData c = new CellData(cellT)
        //                {
        //                    Index = Input_Int32(strs[lut[0]]),
        //                    Location = new Point3D(
        //                        Input_Int32(strs[lut[4]]),
        //                        Input_Int32(strs[lut[5]]),
        //                        Input_Int32(strs[lut[6]])),
        //                };
        //                if (lut[7] > -1) { c.Cell_N = Input_List_Int32(strs[lut[7]]); }
        //                if (lut[8] > -1) { c.Cell_S = (CellData.STATE)Input_Int32(strs[lut[8]]); }
        //                //if (lut[9] > -1) { c.Cell_T = Input_SByte(strs[lut[9]]); }
        //                if (lut[10] > -1) { c.Time_age = Input_Double(strs[lut[10]]); }
        //                if (lut[11] > -1) { c.Time_d = Input_Double(strs[lut[11]]); }
        //                if (lut[12] > -1) { c.Time_m = Input_Double(strs[lut[12]]); }
        //                if (lut[13] > -1) { c.V_m = Input_Double(strs[lut[13]]); }
        //                if (lut[14] > -1) { c.E_AJ = Input_Double(strs[lut[14]]); }
        //                if (lut[15] > -1) { c.E_TJ = Input_Double(strs[lut[15]]); }
        //                if (lut[16] > -1) { c.E_cs = Input_Double(strs[lut[16]]); }
        //                if (lut[17] > -1) { c.E_max = Input_Double(strs[lut[17]]); }
        //                if (lut[18] > -1) { c.Dir_am = (Direction.DIR)Input_Int32(strs[lut[18]]); }
        //                if (lut[19] > -1) { c.Dir_pm = Input_DIR(strs[lut[19]]); }
        //                if (lut[20] > -1) { c.Dir_pd = Input_DIR(strs[lut[20]]); }
        //                if (lut[21] > -1) { c.THETA = Input_Array_Int32(strs, lut[21]); }
        //                if (lut[22] > -1) { c.Time_dev_act = (int)Math.Round(Input_Double(strs[lut[22]]) / CultureTime.Time_step); }
        //                if (lut[23] > -1) { c.Time_dev_total = (int)Math.Round(Input_Double(strs[lut[23]]) / CultureTime.Time_step); }
        //                if (lut[24] > -1 && lut[25] > -1 && lut[26] > -1 && lut[27] > -1 && lut[28] > -1 && lut[29] > -1 && lut[30] > -1) // 2020.03.09
        //                {
        //                    CellMovement.SetValue(
        //                        Input_Double(strs[lut[24]]),
        //                        Input_Double(strs[lut[25]]),
        //                        Input_Double(strs[lut[26]]),
        //                        Input_Double(strs[lut[27]]),
        //                        (int)Math.Round(Input_Double(strs[lut[28]]) / CultureTime.Time_step),
        //                        (int)Math.Round(Input_Double(strs[lut[29]]) / CultureTime.Time_step),
        //                        (int)Math.Round(Input_Double(strs[lut[30]]) / CultureTime.Time_step));
        //                }
        //                while (cells.Count < c.Index)
        //                {
        //                    cells.Add(null);
        //                    CellMovement.SetValue(); // 2020.03.09
        //                }
        //                cells.Add(c);
        //                if (CultureSpace.Xsize - 2 < c.Location.X / 2 && CultureSpace.Ysize - 2 < c.Location.Y)
        //                {
        //                    MessageBox.Show("Out of memory.\nMap size is not correct.");
        //                    //CellData.CellCount = 0;
        //                    return null;
        //                }
        //                CultureSpace.SetMap(c.Location, c.Index);
        //            }
        //            //CellData.CellCount = cells.Count;
        //            return cells;
        //        }
        //        //CellData.CellCount = 0;
        //        return null;
        //    }
        //}

        //private static int[] Input_Array_Int32(string[] strs, int label)
        //{
        //    if (label > 0)
        //    {
        //        string[] s = strs[label].Split('_');
        //        if (s.Length == 20)
        //        {
        //            int[] val = new int[20];
        //            for (int i = 0; i < 20; i++)
        //            {
        //                if (double.TryParse(s[i], out double result))
        //                { val[i] = (int)(result / CultureTime.Time_step); }
        //            }
        //            return val;
        //        }
        //    }
        //    return new int[20];
        //}
        //private static List<Direction.DIR> Input_DIR(string str)
        //{
        //    string[] s = str.Split('_');
        //    List<Direction.DIR> val = new List<Direction.DIR>();
        //    for (int i = 0; i < s.Length; i++)
        //    { if (int.TryParse(s[i], out int result)) { val.Add((Direction.DIR)result); } }
        //    return val;
        //}
        //private static List<int> Input_List_Int32(string str)
        //{
        //    string[] s = str.Split('_');
        //    List<int> val = new List<int>();
        //    for (int i = 0; i < s.Length; i++)
        //    { if (int.TryParse(s[i], out int result)) { val.Add(result); } }
        //    return val;
        //}
        //private static sbyte Input_SByte(string str)
        //{
        //    if (sbyte.TryParse(str, out sbyte result))
        //    { return result; }
        //    return -1;
        //}
        //private static int Input_Int32(string str)
        //{
        //    if (int.TryParse(str, out int result))
        //    { return result; }
        //    return -1;
        //}
        //private static double Input_Double(string str)
        //{
        //    if (double.TryParse(str, out double result))
        //    { return result; }
        //    else if (str == "Inf")
        //    { return double.PositiveInfinity; }
        //    else { return double.NaN; }
        //}
        //#endregion
        #region Random
        private static List<CellData> CellInitialization_Random()
        {
            List<CellData> cells = new List<CellData>();
            int Nsub = GetNsub();
            // 1細胞の平均面積
            double Ac = CultureSpace.Size_lc * CultureSpace.Size_lc * 1E-8 * Math.Sqrt(3) / 2.0; //unit:cm^2
            int ind = 0;
            for (sbyte j = 0; j < X_0.Count; j++)
            {
                int Nseed = (int)(Ac * Nsub * X_0[j]);
                //int n = (int)(AfterSeeding.Get_alpha(j) * Nseed);
                for (int i = 0; i < Nseed; i++)
                {
                    CellData c = new CellData
                    {
                        Index = ind++ // 2021.08.03
                    };
                    cells.Add(c);
                }
                if (Nsub < cells.Count)
                {
                    //CellData.CellCount = 0;
                    return null;
                }
            }
            //CellData.CellCount = cells.Count;
            return cells;
        }
        private static int GetNsub()
        {
            int Nsub = 0;// (CultureSpace.Xsize - 2) * (CultureSpace.Ysize - 2);
            for (int j = 0; j < CultureSpace.Ysize; j++)
            {
                for (int i = 0; i < CultureSpace.Xsize; i++)
                {
                    if (CultureSpace.GetMap(i * 2, j, 1) == -1)
                    {
                        Nsub++;
                    }
                }
            }

            return Nsub;
        }

        private static void InitialCellPlacement_Random(List<CellData> cells)
        {
            for (int i = 0; i < cells.Count; i++)
            {
                while (true)
                {
                    int x = (int)((CultureSpace.Xsize - 2) * Common.Rand_NextDouble()) * 2;
                    int y = (int)((CultureSpace.Ysize - 2) * Common.Rand_NextDouble()) + 1;
                    if (y % 2 == 1)
                    { x++; }
                    Point3D p = new Point3D(x, y, 1);
                    if (CultureSpace.GetMap(p) == -1)
                    {
                        cells[i].Location = p;
                        CultureSpace.SetMap(p, i);
                        break;
                    }
                }
            }
        }
        private static void InitialCellPlacement_Random_CenterBiased(List<CellData> cells) // 2021.08.02
        {
#if false // 正規分布で中心にバイアスをかける
            int x2 = CultureSpace.Xsize / 2;
            double y2 = CultureSpace.Ysize / 2;
            double sigma = Math.Min(x2, y2 * Delta.Cf_y) / 3.0;

            for (int i = 0; i < cells.Count; i++)
            {
                while (true)
                {
                    int x = (int)(x2 + RandomBoxMuller.Next(0, sigma)) * 2;
                    if (x < 0 || CultureSpace.Xsize * 2 - 1 < x) { continue; }
                    int y = (int)(y2 + RandomBoxMuller.Next(0, sigma) / Delta.Cf_y) + 1;
                    if (y < 0 || CultureSpace.Ysize - 1 < y) { continue; }
                    if (y % 2 == 1)
                    { x++; }
                    Point3D p = new Point3D(x, y, 1);
                    if (CultureSpace.GetMap(p) == -1)
                    {
                        cells[i].Location = p;
                        CultureSpace.SetMap(p, i);
                        break;
                    }
                }
            }
#elif false // 円の集中度分布でバイアスをかける
            int y2 = CultureSpace.Ysize / 2 - 1;
            double x2 = y2 % 2 == 0 ? CultureSpace.Xsize / 2 - 1 : CultureSpace.Xsize / 2 - 0.5;
            double y3 = CultureSpace.Ysize * Delta.Cf_y / 2.0;
            double r = (y3 < x2 ? y3 : x2) - 1; // 半径はXとYで小さい方に合わせる

            int[] inds = Common.RandomlySort(cells);
            //for (int i = 0; i < cells.Count; i++)
            foreach (int i in inds)
            {
                while (true)
                {
                    double theta = 2.0 * Math.PI * Common.Rand_NextDouble();
                    // ConcentrationParameter = 0.5 のとき円内に一様分布
                    double r1 = Math.Pow(Common.Rand_NextDouble(), ConcentrationParameter);

                    double _y = r * r1 * Math.Sin(theta);
                    int y = (int)Math.Round(_y / Delta.Cf_y + y2);
                    double _x = r * r1 * Math.Cos(theta);
                    int x = (int)Math.Round(_x + x2) * 2 + (y % 2 == 1 ? 1 : 0);

                    //if (Math.Sqrt(_x * _x + _y * _y) <= r)
                    {
                        Point3D p = new Point3D(x, y, 1);
                        if (CultureSpace.GetMap(p) == -1)
                        {
                            cells[i].Location = p;
                            CultureSpace.SetMap(p, i);
                            break;
                        }
                    }
                    //else { }
                }
            }

#else // 円の集中度分布でバイアスをかける
            // 容器の中心座標
            int y2 = CultureSpace.Ysize / 2 - 1;
            double x2 = y2 % 2 == 0 ? CultureSpace.Xsize / 2 - 1 : CultureSpace.Xsize / 2 - 0.5;
            double y2_ = CultureSpace.Ysize * Delta.Cf_y / 2.0;
            double r = (y2_ < x2 ? y2_ : x2) - 1; // 半径はXとYで小さい方に合わせる

            int[] inds = Common.RandomlySort(cells);
            //for (int i = 0; i < cells.Count; i++)
            foreach (int i in inds)
            {
                while (true)
                {
                    double theta = 2.0 * Math.PI * Common.Rand_NextDouble();// Math.PI / 2.0;//
                    // ConcentrationParameter = 0.5 のとき円内に一様分布
                    double r1 = Math.Pow(Common.Rand_NextDouble(), ConcentrationParameter);

                    double _y = r * r1 * Math.Sin(theta);
                    int y = (int)Math.Round(_y / Delta.Cf_y + y2);

                    //double _x = r * (r1 * Math.Cos(theta) + (r1 - 1.0) * q);
                    double _x = r * (r1 * (Math.Cos(theta) + CenterPeripheryRatio) - CenterPeripheryRatio);
                    int x = (int)Math.Round(_x + x2) * 2 + (y % 2 == 1 ? 1 : 0);

                    {
                        Point3D p = new Point3D(x, y, 1);
                        if (CultureSpace.GetMap(p) == -1)
                        {
                            cells[i].Location = p;
                            CultureSpace.SetMap(p, i);
                            break;
                        }
                    }
                }
            }
#endif
        }
        #endregion
        #region Uniform
        private static List<CellData> CellInitialization_Uniform()
        {
            List<CellData> cells = new List<CellData>();
            int cnt = 0;
            for (int j = 1; j < CultureSpace.Ysize; j += GridInterval)
            {
                for (int i = cnt % 2 == 0 ? 1 : (GridInterval / 2 + 1); i < CultureSpace.Xsize; i += GridInterval)
                {
                    int ii = j % 2 == 0 ? i * 2 : i * 2 + 1;
                    Point3D p = new Point3D(ii, j, 1);
                    if (CultureSpace.GetMap(p) == -1)
                    {
                        int ind = cells.Count;
                        cells.Add(new CellData
                        {
                            Index = ind,
                            Location = p,
                        });
                        CultureSpace.SetMap(p, ind);
                        //CellData.CellCount++;
                    }
                }
                cnt++;
            }
            return cells;
        }
        #endregion
        #region SingleColony_Layer
        private static int SeedingCellNumber_Colony(int layer)
        {
            int y2 = CultureSpace.Ysize / 2 - 1;
            double x2 = y2 % 2 == 0 ? CultureSpace.Xsize / 2 - 1 : CultureSpace.Xsize / 2 - 0.5;

            int cnt = 0;
            for (int j = 1; j < CultureSpace.Ysize - 1; j++)
            {
                for (int i = 1; i < CultureSpace.Xsize - 1; i++)
                {
                    double len = j % 2 == 0
                        ? Delta.GetLength(new Delta((int)((i - x2) * 2), j - y2, 0))
                        : Delta.GetLength(new Delta((int)((i - x2) * 2 + 1), j - y2, 0));
                    if (len <= layer) { cnt++; }
                }
            }
            return cnt;
        }
        private static List<CellData> CellInitialization_SingleColony()
        {
            int n = SeedingCellNumber_Colony(ColonyRadius);
            List<CellData> cells = new List<CellData>();
            //sbyte j = 0;
            for (int i = 0; i < n; i++)
            {
                cells.Add(new CellData());
            }

            //CellData.CellCount = cells.Count;
            return cells;
        }
        private static void InitialCellPlacement_SingleColony(List<CellData> cells)
        {
            //int[] arr = Common.RandomlySort(cells);
            int[] arr = new int[cells.Count];
            for (int i = 0; i < cells.Count; i++)
            { arr[i] = i; }

            int y2 = CultureSpace.Ysize / 2 - 1;
            double x2 = y2 % 2 == 0 ? CultureSpace.Xsize / 2 - 1 : CultureSpace.Xsize / 2 - 0.5;

            int cnt = 0;
            for (int j = 1; j < CultureSpace.Ysize - 1; j++)
            {
                for (int i = 1; i < CultureSpace.Xsize - 1; i++)
                {
                    if (cnt < cells.Count)
                    {
                        double len;
                        int ii;
                        if (j % 2 == 0)
                        {
                            len = Delta.GetLength(new Delta((int)((i - x2) * 2), j - y2, 0));
                            ii = i * 2;
                        }
                        else
                        {
                            len = Delta.GetLength(new Delta((int)((i - x2) * 2 + 1), j - y2, 0));
                            ii = i * 2 + 1;
                        }
                        if (len <= ColonyRadius && CultureSpace.GetMap(i * 2, j, 1) == -1)
                        {
                            int ind = arr[cnt];
                            cells[ind].Index = ind;
                            cells[ind].Location = new Point3D(ii, j, 1);
                            CultureSpace.SetMap(ii, j, 1, ind);
                            cnt++;
                        }
                    }
                }
            }
            if (cnt < cells.Count)
            {
                for (int i = 0; i < cells.Count - cnt; i++)
                {
                    cells.RemoveAt(cells.Count - 1 - i);
                }
            }
        }
        #endregion
        #region SingleColony_CellNumber
        private static List<CellData> CellInitialization_SingleColony_CellNumber()
        {
            List<CellData> cells = new List<CellData>();
            for (int i = 0; i < ColonyCells; i++)
            {
                cells.Add(new CellData());
            }

            //CellData.CellCount = cells.Count;
            return cells;
        }
        private static void InitialCellPlacement_SingleColony_CellNumber(List<CellData> cells)
        {
            int cellNum = (CultureSpace.Xsize - 2) * (CultureSpace.Ysize - 2);
            if (cellNum > cells.Count)
            { cellNum = cells.Count; }

            // 中心座標
            int y = CultureSpace.Ysize / 2;
            int x = y % 2 == 0 ? CultureSpace.Xsize : CultureSpace.Xsize + 1;
            // 中心に最初に配置する // 20191122
            cellNum--;
            Point3D center = BoundaryConditions.Check(x, y, 1);
            cells[0].Index = 0;
            cells[0].Location = center;
            CultureSpace.SetMap(center, 0);


            int k = 1;
            int cnt = 1;
            while (cellNum > 0)
            {
                List<Point3D> Ps = new List<Point3D>();
                for (int j = y - (int)(k / Delta.Cf_y); j <= y + (int)(k / Delta.Cf_y); j++)
                {
                    for (int i = x - k * 2; i <= x + k * 2; i += 2)
                    {
                        // 移動させるy軸の距離によりx軸を補正する
                        int ii = (j - y) % 2 == 0 ? i : i + 1;
                        Point3D p = BoundaryConditions.Check(ii, j, 1);
                        Delta d = new Delta(p.X - x, p.Y - y, 0);
                        if (Delta.GetLength(d) <= k && CultureSpace.GetMap(p) == -1)
                        {
                            Ps.Add(p);
                            cnt++;
                        }
                    }
                }

                if (Ps.Count < cellNum)
                {
                    for (int i = 0; i < Ps.Count; i++)
                    {
                        cellNum--;
                        int ind = cnt - Ps.Count + i;
                        cells[ind].Index = ind;
                        cells[ind].Location = Ps[i];
                        CultureSpace.SetMap(Ps[i], ind);
                    }
                }
                else
                {
                    int tmp = 0;
                    while (cellNum > 0)
                    {
                        int i = (int)(Ps.Count * Common.Rand_NextDouble());
                        if (CultureSpace.GetMap(Ps[i]) == -1)
                        {
                            // 近傍に細胞が存在する場合には配置する
                            for (int j = 7; j < 13; j++) // 20191025
                            {
                                Point3D p = Ps[i] + Delta.GetDelta(j);
                                //Point3D p = BoundaryConditions.Check(Ps[i], Delta.GetDelta(j));
                                if (CultureSpace.GetMap(p) >= 0) // 20191119
                                {
                                    cellNum--;
                                    int ind = cnt - Ps.Count + tmp++;
                                    cells[ind].Index = ind;
                                    cells[ind].Location = Ps[i];
                                    CultureSpace.SetMap(Ps[i], ind);
                                    break;
                                }
                            }
                        }
                    }
                }
                k++;
            }
        }
        #endregion
        #region Multiple colony or aggregation (Random)
        private static List<CellData> CellInitialization_MultipleColony_Random_MeanColony()
        {
            List<CellData> cells = new List<CellData>();
            // 播種コロニー数 // 2020.03.11
            int colonyNum = (int)RandomBoxMuller.Next(NumOfColonies, ColonyNum_SD);
            for (int n = 0; n < colonyNum; n++)
            {
                // 1コロニーに必要な細胞数
                int cellNum = (int)RandomBoxMuller.Next(ColonyCellsMu, ColonyCellsSigma);

                // コロニー中心座標に細胞を播種
                CellInitialization_Multiple_SelectCenter(cells, out int x, out int y, out int z);
                cellNum--; // 播種する残りの細胞数

                int N = 1; // 探索するレイヤー数
                while (cellNum > 0)
                {
                    // 播種する座標の候補を探索
                    List<Point3D> Ps = new List<Point3D>();
                    int k = z;//for (int k = z; k <= z + N; k++)
                    {
                        for (int j = y - (int)(N / Delta.Cf_y); j <= y + (int)(N / Delta.Cf_y); j++)
                        {
                            for (int i = x - N * 2; i <= x + N * 2; i += 2)
                            {
                                // 移動させるy軸の距離によりx軸を補正する
                                int ii = (j - y) % 2 == 0 ? i : i + 1;
                                Point3D p = BoundaryConditions.Check(ii, j, k);
                                // 周期境界条件により移動させる距離が変更されてx軸を補正する必要があるなら
                                if ((j - p.Y) % 2 != 0)
                                { p = BoundaryConditions.Check(ii + 1, p.Y, p.Z); }
                                Delta d = new Delta(ii - x, j - y, k - z);
                                // 範囲内かつ空の単位キューブなら
                                if (Delta.GetLength(d) <= N && CultureSpace.GetMap(p) == -1)
                                { Ps.Add(p); }
                            }
                        }
                    }

                    // 細胞の配置
                    CellInitialization_Multiple_SetCells(cells, ref cellNum, Ps);

                    if (Ps.Count == 0)
                    { N++; }
                }
            }
            //CellData.CellCount = cells.Count;

            return cells;
        }
        private static List<CellData> CellInitialization_MultipleColony_Random_Inoculum()
        {
            List<CellData> cells = new List<CellData>();
            // 基質キューブの数
            int Nsub = GetNsub();
            // 1細胞の平均面積
            double Ac = CultureSpace.Size_lc * CultureSpace.Size_lc * 1E-8 * Math.Sqrt(3) / 2.0; //unit:cm^2
            // 播種する細胞数
            int Nseed = (int)(Ac * Nsub * X_0[0]);
            if (Nseed < NumOfColonies)
            {
                Console.WriteLine("(Nseed:" + Nseed + ") < (NumOfColonies:" + NumOfColonies + ")");
                return null;
            }
            for (int n = 0; n < NumOfColonies; n++)
            {
                // 1コロニーに必要な細胞数
                int cellNum = (int)((double)(n + 1) * Nseed / NumOfColonies) - cells.Count;

                // コロニー中心座標に細胞を播種
                CellInitialization_Multiple_SelectCenter(cells, out int x, out int y, out int z);
                cellNum--; // 播種する残りの細胞数

                int N = 1; // 探索するレイヤー数
                while (cellNum > 0)
                {
                    // 播種する座標の候補を探索
                    List<Point3D> Ps = new List<Point3D>();
                    int k = z;//for (int k = z; k <= z + N; k++)
                    {
                        for (int j = y - (int)(N / Delta.Cf_y); j <= y + (int)(N / Delta.Cf_y); j++)
                        {
                            for (int i = x - N * 2; i <= x + N * 2; i += 2)
                            {
                                // 移動させるy軸の距離によりx軸を補正する
                                int ii = (j - y) % 2 == 0 ? i : i + 1;
                                Point3D p = BoundaryConditions.Check(ii, j, k);
                                // 周期境界条件により移動させる距離が変更されてx軸を補正する必要があるなら
                                if ((j - p.Y) % 2 != 0)
                                { p = BoundaryConditions.Check(ii + 1, p.Y, p.Z); }
                                Delta d = new Delta(ii - x, j - y, k - z);
                                // 範囲内かつ空の単位キューブなら
                                if (Delta.GetLength(d) <= N && CultureSpace.GetMap(p) == -1)
                                { Ps.Add(p); }
                            }
                        }
                    }

                    // 細胞の配置
                    CellInitialization_Multiple_SetCells(cells, ref cellNum, Ps);

                    if (Ps.Count == 0)
                    { N++; }
                }
            }
            //CellData.CellCount = cells.Count;

            return cells;
        }
        private static List<CellData> CellInitialization_MultipleAggregation_Random()
        {
            List<CellData> cells = new List<CellData>();
            for (int n = 0; n < NumOfColonies; n++)
            {
                // 1コロニーに必要な細胞数
                int cellNum = (int)RandomBoxMuller.Next(ColonyCellsMu, ColonyCellsSigma);

                // コロニー中心座標に細胞を播種
                CellInitialization_Multiple_SelectCenter(cells, out int x, out int y, out int z);
                cellNum--; // 播種する残りの細胞数

                int N = 1; // 探索するレイヤー数
                while (cellNum > 0)
                {
                    // 播種する座標の候補を探索
                    List<Point3D> Ps = new List<Point3D>();
                    for (int k = z; k <= z + N; k++)
                    {
                        for (int j = y - (int)(N / Delta.Cf_y); j <= y + (int)(N / Delta.Cf_y); j++)
                        {
                            for (int i = x - N * 2; i <= x + N * 2; i += 2)
                            {
                                // 移動させるy軸の距離によりx軸を補正する
                                int ii = (j - y) % 2 == 0 ? i : i + 1;
                                Point3D p = BoundaryConditions.Check(ii, j, k);
                                // 周期境界条件により移動させる距離が変更されてx軸を補正する必要があるなら
                                if ((j - p.Y) % 2 != 0)
                                { p = BoundaryConditions.Check(ii + 1, p.Y, p.Z); }
                                Delta d = new Delta(ii - x, j - y, k - z);
                                // 範囲内かつ空の単位キューブなら
                                if (Delta.GetLength(d) <= N && CultureSpace.GetMap(p) == -1)
                                { Ps.Add(p); }
                            }
                        }
                    }

                    // 細胞の配置
                    CellInitialization_Multiple_SetCells(cells, ref cellNum, Ps);

                    if (Ps.Count == 0)
                    { N++; }
                }
            }
            //CellData.CellCount = cells.Count;

            return cells;
        }
        private static void CellInitialization_Multiple_SelectCenter(List<CellData> cells, out int x, out int y, out int z)
        {
            z = 1; // 球の中心は基質と接しているものとする
            while (true)
            {
                // 座標
                x = (int)((CultureSpace.Xsize - 2) * Common.Rand_NextDouble()) * 2;
                y = (int)((CultureSpace.Ysize - 2) * Common.Rand_NextDouble()) + 1;
                if (y % 2 == 1) { x++; }

                Point3D p = BoundaryConditions.Check(x, y, z);
                // 選択した座標に細胞が存在していなければ
                if (CultureSpace.GetMap(p) == -1)
                {
                    CellsAdd_MapSet(cells, p);
                    break;
                }
            }
        }
        private static void CellsAdd_MapSet(List<CellData> cells, Point3D p)
        {
            int ind = cells.Count;
            cells.Add(new CellData()
            {
                Index = ind,
                Location = p
            });
            CultureSpace.SetMap(p, ind);
        }
        private static void CellInitialization_Multiple_SetCells(List<CellData> cells, ref int cellNum, List<Point3D> Ps)
        {
            if (Ps.Count < cellNum)
            {
                while (Ps.Count > 0)
                {
                    cellNum--;
                    CellsAdd_MapSet(cells, Ps[0]);
                    Ps.RemoveAt(0);
                }
            }
            else
            {
                while (cellNum > 0 && Ps.Count > 0)
                {
                    int i = (int)(Ps.Count * Common.Rand_NextDouble());
                    if (CultureSpace.GetMap(Ps[i]) == -1)
                    {
                        for (int j = 0; j < 20; j++)
                        {
                            Point3D p = Ps[i] + Delta.GetDelta(j);
                            //Point3D p = BoundaryConditions.Check(Ps[i], Delta.GetDelta(j));
                            // 近傍に細胞が存在する場合には配置する
                            if (CultureSpace.GetMap(p) >= 0)
                            {
                                cellNum--;
                                CellsAdd_MapSet(cells, Ps[i]);
                                Ps.RemoveAt(i);
                                break;
                            }
                        }
                    }
                }
            }
        }
        #endregion
        #region Load cell placement
        //private static List<CellData> CellInitialization_LoadCellPlacement() // 2020.05.27 // 2020.07.18 
        //{
        //    List<CellData> cells = new List<CellData>();
        //    using (Bitmap img = (Bitmap)Image.FromFile(ImagePath))
        //    {
        //        int cnt = 0;
        //        for (int j = 1; j < CultureSpace.Ysize - 1; j++)
        //        {
        //            for (int i = 1; i < CultureSpace.Xsize - 1; i++)
        //            {
        //                int ii = j % 2 == 0 ? i * 2 : i * 2 + 1;
        //                int x = (int)((ii / 2.0 - 0.5) * CultureSpace.Size_lc / PixelSize); // 2020.07.18
        //                int y = (int)((j - 0.5) * Delta.Cf_y * CultureSpace.Size_lc / PixelSize); // 2020.07.18
        //                if (x >= 0 && x < img.Width && y >= 0 && y < img.Height)
        //                {
        //                    int val = img.GetPixel(x, y).R + img.GetPixel(x, y).G + img.GetPixel(x, y).B;
        //                    if (val / 3.0 > 128)
        //                    {
        //                        Point3D p = new Point3D(ii, j, 1);
        //                        if (CultureSpace.GetMap(p) == -1)
        //                        {
        //                            CellData c = new CellData()
        //                            {
        //                                Location = p,
        //                                Index = cnt++
        //                            };
        //                            cells.Add(c);
        //                            CultureSpace.SetMap(p, c.Index);
        //                        }
        //                    }
        //                }
        //            }
        //        }
        //    }
        //    return cells;
        //}
        #endregion Load cell placement
        #region Non-uniform distribution
        private static void InitialCellPlacement_NonUniform(ref List<CellData> cells) // 2021.07.07
        {
            for (int i = 0; i < cells.Count; i++)
            {
                int cnt = 0;
                while (true)
                {
                    int x = (int)((CultureSpace.Xsize - 2) * Common.Rand_NextDouble()) * 2;
                    int y = (int)((CultureSpace.Ysize - 2) * Common.Rand_NextDouble()) + 1;
                    int xmax = (CultureSpace.Xsize - 2) * 2; // 2021.07.09
                    if (y % 2 == 1)
                    { x++; xmax++; }
                    Point3D p = new Point3D(x, y, 1);
                    if (CultureSpace.GetMap(p) == -1)
                    {
                        double px;
                        if (Uniformity_Inf)
                        { px = 1.0; }
                        else
                        {
                            double axmax = Uniformity * xmax; // 2021.07.09
                            px = Math.Exp(-x / axmax) / (axmax * (1.0 - Math.Exp(-1.0 / Uniformity)));
                            if (double.IsNaN(px))
                            {

                            }
                        }
                        if (Common.Rand_NextDouble() < px)
                        {
                            cells[i].Location = p;
                            CultureSpace.SetMap(p, i);
                            break;
                        }
                    }
                    if (cnt++ == int.MaxValue)
                    {
                        cells = null;
                        break;
                    }
                }
                if (cells == null)
                { break; }
            }
        }
        #endregion

    }
    public static class RandomBoxMuller
    {
        /// <summary>
        /// 正規分布(μ・σ)に従う乱数
        /// </summary>
        /// <param name="mu">平均:μ</param>
        /// <param name="sigma">標準偏差:σ</param>
        /// <param name="getCos">Cos or Sin:True->cos</param>
        /// <returns>N(mu,sigma)に従う乱数</returns>
        public static double Next(double mu = 0.0, double sigma = 1.0, bool getCos = true)
        {
            if (getCos)
            {
                double rand;
                while ((rand = Common.Rand_NextDouble()) == 0.0) ;
                double rand2 = Common.Rand_NextDouble();
                double normrand = Math.Sqrt(-2.0 * Math.Log(rand)) * Math.Cos(2.0 * Math.PI * rand2);
                normrand = normrand * sigma + mu;
                return normrand;
            }
            else
            {
                double rand;
                while ((rand = Common.Rand_NextDouble()) == 0.0) ;
                double rand2 = Common.Rand_NextDouble();
                double normrand = Math.Sqrt(-2.0 * Math.Log(rand)) * Math.Sin(2.0 * Math.PI * rand2);
                normrand = normrand * sigma + mu;
                return normrand;
            }
        }

        /// <summary>
        /// 正規分布(μ・σ)に従う乱数
        /// </summary>
        /// <param name="mu">平均:μ</param>
        /// <param name="sigma">標準偏差:σ</param>
        /// <returns>N(mu,sigma)に従う乱数</returns>
        public static double[] NextPair(double mu = 0.0, double sigma = 1.0)
        {
            double[] normrand = new double[2];
            double rand;
            while ((rand = Common.Rand_NextDouble()) == 0.0) ;
            double rand2 = Common.Rand_NextDouble();
            normrand[0] = Math.Sqrt(-2.0 * Math.Log(rand)) * Math.Cos(2.0 * Math.PI * rand2);
            normrand[0] = normrand[0] * sigma + mu;
            normrand[1] = Math.Sqrt(-2.0 * Math.Log(rand)) * Math.Sin(2.0 * Math.PI * rand2);
            normrand[1] = normrand[1] * sigma + mu;
            return normrand;
        }

        /// <summary>
        /// 一様乱数を取得する（暗号化サービス プロバイダーを使用）
        /// </summary>
        /// <returns>乱数:-1～+1</returns>
        public static double GetRandom()
        {
            using (RNGCryptoServiceProvider rng = new RNGCryptoServiceProvider())
            {
                byte[] bs = new byte[sizeof(int)];
                rng.GetBytes(bs);
                int iR = BitConverter.ToInt32(bs, 0);
                return (double)iR / int.MaxValue;
            }
        }

        /// <summary>
        /// 標準正規分布に従う乱数
        /// </summary>
        /// <returns>N(0,1)に従う乱数</returns>
        public static double GetNormRandom()
        {
            double dR1 = Math.Abs(GetRandom());
            double dR2 = Math.Abs(GetRandom());
            return Math.Sqrt(-2 * Math.Log(dR1, Math.E)) * Math.Cos(2 * Math.PI * dR2);
        }

        /// <summary>
        /// 正規分布(μ・σ)に従う乱数
        /// </summary>
        /// <param name="m">平均:μ</param>
        /// <param name="s">標準偏差:σ</param>
        /// <returns>N(m,s)に従う乱数</returns>
        public static double GetNormRandom(double m, double s)
        {
            return m + s * GetNormRandom();
        }
    }
    public static class Output
    {
        private static string OutputDir = "";

        public static void Run(List<CellData> cells, int step = 0)
        {
            //// ゼロなら出力しない
            //if (StepInterval != 0)
            //{
            //if (step % StepInterval == 0)
            //{
            Run_CellData(cells, step);
            //}
            //}
        }

        #region // Run_CellData
        public const string Name_Index = "Index";
        public const string Name_X = "X (mm)";
        public const string Name_Y = "Y (mm)";
        public const string Name_hexX = "hexX (-)";
        public const string Name_hexY = "hexY (-)";
        private static void Run_CellData(List<CellData> cells, int step)
        {
            string name = "cell.csv"; // 2021.02.16
            string path = Path.Combine(OutputDir, string.Format(name, step)); // 2021.02.16
            try
            {
                using StreamWriter sw = new StreamWriter(File.Open(path, FileMode.Create, FileAccess.Write));
                string s = Name_Index;
                s += "," + Name_X + "," + Name_Y + "," + Name_hexX + "," + Name_hexY;
                sw.WriteLine(s);
                for (int i = 0; i < cells.Count; i++)
                {
                    // 死細胞は保存しない
                    if (cells[i] != null)
                    {
                        string str = cells[i].Index.ToString();
                        str += "," + cells[i].Location.X / 2.0 * CultureSpace.Size_lc * 1E-3;
                        str += "," + cells[i].Location.Y * Delta.Cf_y * CultureSpace.Size_lc * 1E-3;
                        str += "," + cells[i].Location.X;
                        str += "," + cells[i].Location.Y;
                        sw.WriteLine(str);
                    }
                }
            }
            catch (Exception ex) // 2020.12.02
            {
                Console.Write(ex.Message + "\r\nRetry? (y/n) ");
                if (Console.Read() == 'y')
                {
                    Run_CellData(cells, step);
                }
            }
        }
        //private static int Digit(int num)
        //{
        //    // Math.Log10(0)はNegativeInfinityを返すため、別途処理する。
        //    return (num == 0) ? 1 : ((int)Math.Log10(num) + 1);
        //}
        //private static string Str_(double val)
        //{
        //    return double.IsInfinity(val) ? "Inf" : val.ToString();
        //}
        //private static string Str_(List<int> Cell_N)
        //{
        //    string str = Cell_N[0].ToString();
        //    for (int j = 1; j < Cell_N.Count; j++)
        //    { str += "_" + Cell_N[j]; }
        //    return str;
        //}
        //private static string Str_(List<Direction.DIR> i)
        //{
        //    if (i.Count > 0)
        //    {
        //        string str = ((int)i[0]).ToString();
        //        for (int j = 1; j < i.Count; j++)
        //        { str += "_" + ((int)i[j]).ToString(); }
        //        return str;
        //    }
        //    else
        //    { return ""; }
        //}
        //private static string Str_(int[] i, double val)
        //{
        //    string str = (i[0] * val).ToString();
        //    for (int j = 1; j < i.Length; j++)
        //    { str += "_" + i[j] * val; }
        //    return str;
        //}
        #endregion
    }
    public class Capturing
    {
        public void Run(double interval)
        {
            double len_u = 2E3; // 1辺の長さ (um)
            int len_p = 1000; // 1辺の長さ (pixel)

            Point3D[] position = new Point3D[25];
            double[,] points = new double[25, 2];
            // 容器の中心座標
            int y2 = CultureSpace.Ysize / 2 - 1;
            double x2 = y2 % 2 == 0 ? CultureSpace.Xsize / 2 - 1 : CultureSpace.Xsize / 2 - 0.5;
            int cnt = 0;
            for (int j = -2; j <= 2; j++)
            {
                for (int i = -2; i <= 2; i++)
                {
                    double _y = (j * interval * 1E3 - len_p) / CultureSpace.Size_lc;
                    int y = (int)Math.Round(_y / Delta.Cf_y + y2);
                    double _x = (i * interval * 1E3 - len_p) / CultureSpace.Size_lc;
                    int x = (int)Math.Round(_x + x2) * 2 + (y % 2 == 1 ? 1 : 0);
                    points[cnt, 0] = _x;
                    points[cnt, 1] = _y;
                    position[cnt++] = new Point3D(x, y, 1);
                }
            }

            int height = (int)Math.Round(len_u / CultureSpace.Size_lc / Delta.Cf_y);
            int width = (int)Math.Round(len_u / CultureSpace.Size_lc);

            for (int k = 0; k < 25; k++)
            {
                Bitmap bmp = new(len_p, len_p);
                for (int j = position[k].Y; j < position[k].Y + height; j++)
                {
                    for (int i = position[k].X / 2; i < position[k].X / 2 + width; i++)
                    {
                        int val = CultureSpace.GetMap(i * 2, j, 1);
                        if (val >= 0)
                        {

                        }
                    }
                }
            }
        }
    }
}
