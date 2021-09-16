using System;
using System.Collections.Generic;
using System.IO;
using System.Text;
using System.Drawing;

namespace SpatialHeterogeneity
{
    class Program
    {
        static void Main(string[] args)
        {
#if false
            double A_c = 1.4e-6; // Average cell area (cm^2)
            lc = Math.Sqrt(2.0 / Math.Sqrt(3.0) * A_c) * 1E4; // (um)

            radius = 11 * 1E3; // Radius of culture vessel (um)
            double v = radius * 2 / lc;
            xSize = (int)(v + 2);
            ySize = (int)(v / Delta.Cf_y + 2);

            density = 5e3; // Seeding density (cells/cm^2)

            interval = 1.9 * 1E3; // Image capture interval (um)

            N_col = new int[1] { 1 };
            m_con = new double[1] { 0.5 };
            p_con = new double[1] { 0.0 };

            //output = @"C:\Users\foo\Desktop";
            output = AppDomain.CurrentDomain.BaseDirectory;
#else
            string openfilename = "";
            if (args.Length > 0)
            {
                openfilename = args[0];
            }
            else if (File.Exists("Parameter.csv"))
            {
                openfilename = "Parameter.csv";
            }
            if (openfilename == "" || !ReadFile(openfilename))
            {
                return;
            }
#endif
            // progress
            Console.CursorVisible = false;
            char[] bars = { '／', '―', '＼', '｜' };
            int sum = N_col.Length * m_con.Length * p_con.Length;
            int cnt = 0;

            for (int i = 0; i < N_col.Length; i++)
            {
                for (int j = 0; j < m_con.Length; j++)
                {
                    for (int k = 0; k < p_con.Length; k++)
                    {
                        // 回転する棒を表示
                        Console.Write(bars[cnt++ % 4]);
                        // 進むパーセンテージを表示
                        Console.Write("{0, 4:d0}% ({1} / {2})", 100 * cnt / sum, cnt, sum);
                        // カーソル位置を初期化
                        Console.SetCursorPosition(0, Console.CursorTop);

                        BoundaryConditions.SetParameter(false, false, false);
                        CultureSpace.SetParameter(xSize, ySize, 3, lc, 1);
                        // マップ作成
                        if (!CultureSpace.MapCreation())
                        { return; }
                        // 細胞の配置
                        if (!Seeding.Run_Random_CenterBiased_Colony(
                            out List<CellData> cells, density, m_con[j], p_con[k], N_col[i]))
                        { return; }

                        string folder = "{N_col=" + N_col[i] + "}{m_con=" + m_con[j] + "}{p_con=" + p_con[k] + "}";
                        Capturing capturing = new(Path.Combine(output, folder), interval, radius);
                        capturing.Run();

                        // 保存
                        //Output.Run(cells);
                    }
                }
            }

            Console.CursorVisible = true;
        }
        static double lc; // (um)
        static double radius; // (um)
        static int xSize;
        static int ySize;
        static double density; // (cells/cm^2)
        static double interval; // (um)
        static int[] N_col;
        static double[] m_con;
        static double[] p_con;
        static string output;

        private static bool ReadFile(string filename)
        {
            Encoding.RegisterProvider(CodePagesEncodingProvider.Instance);
            var sjis = Encoding.GetEncoding("shift_jis");
            using (FileStream fs = new(filename, FileMode.Open, FileAccess.Read, FileShare.ReadWrite))
            using (StreamReader sr = new(fs, sjis))
            {
                try
                {
                    if (ReadLine(sr, "A_c", out double Ac))
                    {
                        lc = Math.Sqrt(2.0 / Math.Sqrt(3.0) * Ac) * 1E4;
                        //Console.WriteLine("Distance between cells l_c = " + lc + " (um)");
                    }
                    else return false;
                    if (ReadLine(sr, "R_vessel", out double r))
                    {
                        radius = r * 1E3; // mm -> um
                        double v = radius * 2 / lc;
                        xSize = (int)(v + 2);
                        ySize = (int)(v / Delta.Cf_y + 2);
                        //Console.WriteLine("Xsize = " + xSize + " (-)");
                        //Console.WriteLine("Ysize = " + ySize + " (-)");
                    }
                    else return false;
                    if (!ReadLine(sr, "X_0", out density)) return false;
                    if (!ReadLine(sr, "L_interval", out interval)) return false;
                    interval *= 1E3; // (um)
                    sr.ReadLine();
                    if (!ReadLine(sr, "N_colony", out N_col)) return false;
                    if (!ReadLine(sr, "m_con", out m_con)) return false;
                    if (!ReadLine(sr, "p_con", out p_con)) return false;
                    string[] s = sr.ReadLine().Split(',');
                    if (s[0] == "outputdirectory") output = s[1];
                    else output = AppDomain.CurrentDomain.BaseDirectory;
                    if (!Directory.Exists(output))
                    {
                        Console.WriteLine("OutputDirectory: " + AppDomain.CurrentDomain.BaseDirectory);
                        output = AppDomain.CurrentDomain.BaseDirectory;
                    }
                }
                catch (Exception ex)
                {
                    Console.WriteLine(ex.Message);
                    Console.WriteLine("Exception in ReadFile");
                    return false;
                }
            }
            return true;
        }
        private static bool ReadLine(StreamReader sr, string str, out double val)
        {
            string[] strs = sr.ReadLine().Split(',');
            if (strs[0] == str)
            {
                return double.TryParse(strs[1], out val);
            }
            else
            {
                val = double.NaN;
                return false;
            }
        }
        private static bool ReadLine(StreamReader sr, string str, out int val)
        {
            string[] strs = sr.ReadLine().Split(',');
            if (strs[0] == str)
            {
                return int.TryParse(strs[1], out val);
            }
            else
            {
                val = 0;
                return false;
            }
        }
        private static bool ReadLine(StreamReader sr, string str, out double[] val)
        {
            string[] strs = sr.ReadLine().Split(',');
            if (strs[0] == str)
            {
                if (double.TryParse(strs[1], out double vmin) &&
                    double.TryParse(strs[2], out double vmax) &&
                    decimal.TryParse(strs[3], out decimal vinc))
                {
                    int precision = GetPrecision(vinc);
                    double dinc = (double)vinc;
                    double a = Math.Pow(10, precision);
                    vmin *= a;
                    vmax *= a;
                    dinc *= a;
                    int min = (int)Math.Round(vmin);
                    int max = (int)Math.Round(vmax);
                    int inc = (int)Math.Round(dinc);
                    int cnt = vinc == 0 ? 1 : (max - min) / inc + 1;
                    val = new double[cnt];
                    val[0] = min / a;
                    for (int i = 1; i < cnt; i++)
                    {
                        val[i] = (min + inc * i) / a;
                    }
                    return true;
                }
            }
            val = null;
            return false;
        }
        private static bool ReadLine(StreamReader sr, string str, out int[] val)
        {
            string[] strs = sr.ReadLine().Split(',');
            if (strs[0] == str)
            {
                if (int.TryParse(strs[1], out int vmin) &&
                    int.TryParse(strs[2], out int vmax) &&
                    int.TryParse(strs[3], out int vinc))
                {
                    int cnt = vinc == 0 ? 1 : (vmax - vmin) / vinc + 1;
                    val = new int[cnt];
                    val[0] = vmin;
                    for (int i = 1; i < cnt; i++)
                    {
                        val[i] = val[i - 1] + vinc;
                    }
                    return true;
                }
            }
            val = null;
            return false;
        }
        /// <summary>
        /// 小数点以下の桁数を取得
        /// </summary>
        private static int GetPrecision(decimal price)
        {
            string priceString = price.ToString().TrimEnd('0');

            int index = priceString.IndexOf('.');
            if (index == -1)
                return 0;

            return priceString.Substring(index + 1).Length;
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
            Random rand = new();
            for (int i = 0; i < arr.Count; i++)
            {
                int val = (int)(rand.NextDouble() * (arr.Count - i)) + i;
                int tmp = arr[i];
                arr[i] = arr[val];
                arr[val] = tmp;
            }
            return arr.ToArray();
        }
    }
    public class CellData
    {
        /// <summary>
        /// CellData の新規作成
        /// </summary>
        public CellData() // 2020.12.03
        {
            // 細胞インデックス
            Index = -1;
            // 細胞位置
            Location = new Point3D(-1, -1, -1);
        }

        /// <summary>
        /// 細胞のインデックスを取得または設定します。
        /// </summary>
        internal int Index { get; set; }
        /// <summary>
        /// 細胞の座標を取得または設定します。
        /// </summary>
        internal Point3D Location { get; set; }
    }
    public class Point3D
    {
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

        public static Point3D operator +(Point3D left, Delta right)
        {
            return BoundaryConditions.Check(left, right);
        }
    }
    public class Delta : Point3D
    {
        public Delta(int dx, int dy, int dz) : base(dx, dy, dz)
        { DX = dx; DY = dy; DZ = dz; }

        public int DX { get; set; }
        public int DY { get; set; }
        public int DZ { get; set; }
        /// <summary>
        /// 六角格子座標から実座標へのy軸の補正係数Correction factor: sqrt(3)/2
        /// </summary>
        public static readonly double Cf_y = Math.Sqrt(3) / 2.0;

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
        private static TYPE _type;
        private static int[,,] Map;

        // 静的プロパティ
        internal static double Size_lc { get; private set; }
        internal static int Xsize { get; private set; }
        internal static int Ysize { get; private set; }
        internal static int Zsize { get; private set; }

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
        private static void Creation_Circle() // 2021.09.15
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
                        for (int k = 1; k < Zsize - 1; k++)
                        {
                            SetMap(i * 2, j, k, -1); // empty grid
                        }
                    }
                    else
                    {
                        for (int k = 1; k < Zsize - 1; k++)
                        {
                            SetMap(i * 2, j, k, -2); // block
                        }
                    }
                    SetMap(i * 2, j, 0, -3); // substrate
                    SetMap(i * 2, j, Zsize - 1, -2); // block
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

        internal static void SetMap(int row, int col, int dep, int value)
        {
            Map[dep, col, row / 2] = value;
        }
        internal static void SetMap(Point3D point, int value)
        {
            SetMap(point.X, point.Y, point.Z, value);
        }
    }
    public class Seeding
    {
        private static double X_0; // Inoculum size 接種細胞密度
        private static int ColonyCells; // コロニー内の細胞数
        private static double m_con; // 2021.08.02
        private static double p_con; // 2021.08.31

        // 実処理
        public static bool Run_Random_CenterBiased_Colony(out List<CellData> cells, double X0, int colonyCells = 1)
        {
            X_0 = X0;
            m_con = 0.5;
            p_con = 0.0;
            ColonyCells = colonyCells;

            cells = CellInitialization_Random();
            if (cells == null) { return false; }
            InitialCellPlacement_Random_CenterBiased_Colony(cells);
            return cells != null;
        }
        public static bool Run_Random_CenterBiased_Colony(out List<CellData> cells, double X0, double concentrationParameter = 0.5, double centerPeripheryRatio = 0.0, int colonyCells = 1)
        {
            X_0 = X0;
            m_con = concentrationParameter;
            p_con = centerPeripheryRatio;
            ColonyCells = colonyCells;

            cells = CellInitialization_Random();
            if (cells == null) { return false; }
            InitialCellPlacement_Random_CenterBiased_Colony(cells);
            return cells != null;
        }

        #region Random
        private static List<CellData> CellInitialization_Random()
        {
            List<CellData> cells = new List<CellData>();
            int Nsub = GetNsub();
            // 1細胞の平均面積
            double Ac = CultureSpace.Size_lc * CultureSpace.Size_lc * 1E-8 * Math.Sqrt(3) / 2.0; //unit:cm^2
            int ind = 0;
            int Nseed = (int)(Ac * Nsub * X_0);
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
                return null;
            }
            return cells;
        }
        private static int GetNsub()
        {
            int Nsub = 0;
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

        private static void InitialCellPlacement_Random_CenterBiased_Colony(List<CellData> cells) // 2021.08.02
        {
            // 円の集中度分布でバイアスをかける

            // 容器の中心座標
            int y2 = CultureSpace.Ysize / 2 - 1;
            double x2 = y2 % 2 == 0 ? CultureSpace.Xsize / 2 - 1 : CultureSpace.Xsize / 2 - 0.5;
            double y2_ = CultureSpace.Ysize * Delta.Cf_y / 2.0;
            double r = (y2_ < x2 ? y2_ : x2) - 1; // 半径はXとYで小さい方に合わせる

            int[] inds = Common.RandomlySort(cells);

            Random rand = new();

            for (int ind = 0; ind < cells.Count; ind++)
            {
                while (true)
                {
                    double theta = 2.0 * Math.PI * rand.NextDouble();// Math.PI / 2.0;//
                    // m_con = 0.5 のとき円内に一様分布
                    double r1 = Math.Pow(rand.NextDouble(), m_con);

                    double _y = r * r1 * Math.Sin(theta);
                    int y = (int)Math.Round(_y / Delta.Cf_y + y2);

                    //double _x = r * (r1 * Math.Cos(theta) + (r1 - 1.0) * q);
                    double _x = r * (r1 * (Math.Cos(theta) + p_con) - p_con);
                    int x = (int)Math.Round(_x + x2) * 2 + (y % 2 == 1 ? 1 : 0);

                    Point3D p = new Point3D(x, y, 1);
                    if (CultureSpace.GetMap(p) == -1)
                    {
                        // 1コロニーに必要な細胞数
                        int cellNum = ColonyCells;

                        // 中心座標に細胞を配置
                        cellNum--;
                        if (cells.Count <= ++ind) break;
                        cells[inds[ind]].Location = p;
                        CultureSpace.SetMap(p, inds[ind]);

                        int N = 1; // 探索するレイヤー数
                        while (cellNum > 0 && ind < cells.Count)
                        {
                            // 播種する座標の候補を探索
                            List<Point3D> Ps = new List<Point3D>();
                            for (int j = y - (int)(N / Delta.Cf_y); j <= y + (int)(N / Delta.Cf_y); j++)
                            {
                                for (int i = x - N * 2; i <= x + N * 2; i += 2)
                                {
                                    // 移動させるy軸の距離によりx軸を補正する
                                    int ii = (j - y) % 2 == 0 ? i : i + 1;
                                    Point3D pp = BoundaryConditions.Check(ii, j, 1);
                                    // 周期境界条件により移動させる距離が変更されてx軸を補正する必要があるなら
                                    if ((j - pp.Y) % 2 != 0)
                                    { pp = BoundaryConditions.Check(ii + 1, pp.Y, pp.Z); }
                                    Delta d = new Delta(ii - x, j - y, 0);
                                    // 範囲内かつ空の単位キューブなら
                                    if (Delta.GetLength(d) <= N && CultureSpace.GetMap(pp) == -1)
                                    { Ps.Add(pp); }
                                }
                            }

                            // 細胞の配置
                            if (Ps.Count < cellNum)
                            {
                                while (Ps.Count > 0)
                                {
                                    cellNum--;
                                    if (cells.Count <= ++ind) break;
                                    cells[inds[ind]].Location = Ps[0];
                                    CultureSpace.SetMap(Ps[0], inds[ind]);
                                    Ps.RemoveAt(0);
                                }
                            }
                            else
                            {
                                while (cellNum > 0 && Ps.Count > 0 && ind < cells.Count)
                                {
                                    int i = (int)(Ps.Count * rand.NextDouble());
                                    if (CultureSpace.GetMap(Ps[i]) == -1)
                                    {
                                        for (int j = 7; j < 13; j++)
                                        {
                                            Point3D pp = Ps[i] + Delta.GetDelta(j);
                                            // 近傍に細胞が存在する場合には配置する
                                            if (CultureSpace.GetMap(pp) >= 0)
                                            {
                                                cellNum--;
                                                if (cells.Count <= ++ind) break;
                                                cells[inds[ind]].Location = Ps[i];
                                                CultureSpace.SetMap(Ps[i], inds[ind]);
                                                Ps.RemoveAt(i);
                                                break;
                                            }
                                        }
                                    }
                                }
                            }

                            if (Ps.Count == 0)
                            { N++; }
                        }

                        break;
                    }
                }
            }
        }
        #endregion
    }
    public static class Output
    {
        private static string OutputDir = "";

        public static void Run(List<CellData> cells, int step = 0)
        {
            Run_CellData(cells, step);
        }

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
    }
    public class Capturing
    {
        public Capturing(string path, double interval, double radius)
        {
            if (!Directory.Exists(path))
            {
                Directory.CreateDirectory(path);
            }
            OutputDirectory = path;
            Interval = interval;
            Radius = radius;
        }
        private string OutputDirectory;
        private double Interval; // (um)
        private double Radius; // (um)

        /// <summary>
        /// 
        /// </summary>
        /// <param name="len_um">1辺の長さ (um)</param>
        /// <param name="len_px">1辺の長さ (pixel)</param>
        public void Run(double len_um = 2E3, int len_px = 1000)
        {
            double res = len_px / len_um; // 解像度

            // 容器の中心座標
            int y_center = CultureSpace.Ysize / 2 - 1;
            double x_center = y_center % 2 == 0 ? CultureSpace.Xsize / 2 - 1 : CultureSpace.Xsize / 2 - 0.5;
            int cnt = 0;
            // 25点の位置を計算
            Point[] position = new Point[25]; // 左上の座標
            for (int j = -2; j <= 2; j++)
            {
                for (int i = -2; i <= 2; i++)
                {
                    double _y = j * Interval - len_um / 2; // (um)
                    int y = (int)Math.Floor(_y / CultureSpace.Size_lc / Delta.Cf_y + y_center); // (-)
                    if (y % 2 == 1) { y++; }
                    double _x = i * Interval - len_um / 2; // (um)
                    int x = (int)Math.Floor(_x / CultureSpace.Size_lc + x_center) * 2 + (y % 2 == 1 ? 1 : 0); // (-)
                    position[cnt++] = new Point(x, y);
                }
            }

            int height = (int)Math.Round(len_um / CultureSpace.Size_lc / Delta.Cf_y); // (-)
            int width = (int)Math.Round(len_um / CultureSpace.Size_lc); // (-)

            for (int k = 0; k < 25; k++)
            {
                Bitmap bmp = new(len_px, len_px);
                Graphics graphics = Graphics.FromImage(bmp);
                graphics.FillRectangle(Brushes.Black, graphics.VisibleClipBounds);
                graphics.Dispose();
                for (int j = position[k].Y; j < position[k].Y + height; j++)
                {
                    for (int i = position[k].X; i < position[k].X + width * 2; i += 2)
                    {
                        int val = CultureSpace.GetMap(i, j, 1);
                        if (val >= 0)
                        {
                            int x = (int)((i - position[k].X) / 2.0 * CultureSpace.Size_lc * res); // (px)
                            int y = (int)((j - position[k].Y) * Delta.Cf_y * CultureSpace.Size_lc * res); // (px)
                            if (x >= 0 && x < len_px && y >= 0 && y < len_px)
                            {
                                bmp.SetPixel(x, y, Color.White);
                            }
                        }
                    }
                }
                bmp.Save(Path.Combine(OutputDirectory, string.Format("{0:00}", k + 1) + ".bmp"));
            }

            // 容器全体
            {
                int mx = 3;
                Bitmap img = new(CultureSpace.Xsize * mx, (int)(CultureSpace.Ysize * Delta.Cf_y * mx));
                Graphics g = Graphics.FromImage(img);
                // 全体をグレーにする
                g.FillRectangle(Brushes.Gray, g.VisibleClipBounds);
                // 容器内を黒で描画
                float hf = (float)(Radius / CultureSpace.Size_lc * mx);
                float xf = img.Size.Width / 2.0f - hf;
                float yf = img.Size.Height / 2.0f - hf;
                float wf = hf * 2f;
                g.FillEllipse(Brushes.Black, xf, yf, wf, wf);
                // 細胞を描画
                for (int j = 0; j < CultureSpace.Ysize; j++)
                {
                    for (int i = 0; i < CultureSpace.Xsize; i++)
                    {
                        int ii = j % 2 == 0 ? i * 2 : i * 2 + 1;
                        if (CultureSpace.GetMap(new Point3D(ii, j, 1)) >= 0)
                        {
                            int x = ii * mx / 2;
                            int y = (int)(j * mx * Delta.Cf_y);
                            g.FillEllipse(Brushes.White, new RectangleF(x - mx / 2, y - mx / 2, mx, mx));
                            g.DrawEllipse(Pens.Gray, x - mx / 2, y - mx / 2, mx, mx);
                        }
                    }
                }
                // 赤枠を付ける
                for (int k = 0; k < 25; k++)
                {
                    int j = position[k].Y;
                    int i = position[k].X;
                    int ii = j % 2 == 0 ? i * 2 : i * 2 + 1;
                    int x = ii * mx / 2 / 2;
                    int y = (int)(j * mx * Delta.Cf_y);
                    g.DrawRectangle(Pens.Red, x, y, width * mx, (float)(height * mx * Delta.Cf_y));
                }
                g.Dispose();
                img.Save(Path.Combine(OutputDirectory, "vessel.bmp"));
            }

        }
    }
}
