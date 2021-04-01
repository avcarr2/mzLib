﻿// Copyright 2016 Stefan Solntsev
//
// This file (PeriodicTable.cs) is part of UsefulProteomicsDatabases.
//
// UsefulProteomicsDatabases is free software: you can redistribute it and/or modify it
// under the terms of the GNU Lesser General Public License as published
// by the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// UsefulProteomicsDatabases is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
// License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with UsefulProteomicsDatabases. If not, see <http://www.gnu.org/licenses/>.

using Chemistry;
using MzLibUtil;
using System;
using System.Globalization;
using System.Text.RegularExpressions;

namespace UsefulProteomicsDatabases
{
    /// <summary>
    /// The Periodic Table of Elements.
    /// </summary>
    public static class PeriodicTableLoader
    {
        public static void Load()
        {
            var periodicTable = thePeriodicTable.Split(new string[] { "\n" }, StringSplitOptions.None);

            int? atomicNumber = null;
            string atomicSymbol = null;
            int? isotopeNumber = null;
            double? atomicMass = null;
            double? abundance = null;
            double? averageMass = null;

            foreach (string line in periodicTable)
            {
                if (string.IsNullOrWhiteSpace(line))
                {
                    continue;
                }

                string key = line.Substring(0, line.IndexOf('=')).Trim();
                string value = line.Substring(line.IndexOf('=') + 1).Trim();

                switch (key)
                {
                    case "Atomic Number":
                        atomicNumber = int.Parse(value);
                        break;

                    case "Atomic Symbol":
                        atomicSymbol = value;
                        break;

                    case "Mass Number":
                        isotopeNumber = int.Parse(value);
                        break;

                    case "Relative Atomic Mass":
                        atomicMass = double.Parse(Regex.Match(line, @"[\d\.]+").Value, CultureInfo.InvariantCulture);
                        break;

                    case "Isotopic Composition":
                        if (Regex.Match(line, @"[\d\.]+").Success)
                        {
                            abundance = double.Parse(Regex.Match(line, @"[\d\.]+").Value, CultureInfo.InvariantCulture);
                        }
                        break;

                    case "Standard Atomic Weight":
                        if (Regex.Match(line, @"\[").Success)
                        {
                            double averageMass1 = double.Parse(Regex.Match(line, @"(?<=\[)[\d\.]+").Value, CultureInfo.InvariantCulture);
                            string averageMass2StringRep = Regex.Match(line, @"(?<=,)[\d\.]+").Value;

                            if (double.TryParse(averageMass2StringRep, out double averageMass2))
                            {
                                averageMass = (averageMass1 + averageMass2) / 2;
                            }
                            else
                            {
                                averageMass = averageMass1;
                            }
                        }
                        else if (Regex.Match(line, @"[\d\.]+").Success)
                        {
                            averageMass = Convert.ToDouble(Regex.Match(line, @"[\d\.]+").Value, CultureInfo.InvariantCulture);
                        }
                        break;

                    case "Notes":
                        if (PeriodicTable.GetElement(atomicNumber.Value) == null)
                        {
                            if (atomicSymbol != null && atomicNumber.HasValue && averageMass.HasValue && isotopeNumber.HasValue && atomicMass.HasValue && abundance.HasValue)
                            {
                                var element = new Element(atomicSymbol, atomicNumber.Value, averageMass.Value);
                                element.AddIsotope(isotopeNumber.Value, atomicMass.Value, abundance.Value);
                                PeriodicTable.Add(element);
                            }
                        }
                        else if (abundance.HasValue)
                        {
                            PeriodicTable.GetElement(atomicNumber.Value).AddIsotope(isotopeNumber.Value, atomicMass.Value, abundance.Value);
                        }

                        atomicNumber = null;
                        atomicSymbol = null;
                        isotopeNumber = null;
                        atomicMass = null;
                        abundance = null;
                        averageMass = null;
                        break;

                    default:
                        throw new MzLibException("Could not parse line from periodic table: " + line);
                }
            }
        }

        private const string thePeriodicTable =
            @"
            Atomic Number = 1
            Atomic Symbol = H
            Mass Number = 1
            Relative Atomic Mass = 1.00782503223(9)
            Isotopic Composition = 0.999885(70)
            Standard Atomic Weight = [1.00784,1.00811]
            Notes = m

            Atomic Number = 1
            Atomic Symbol = D
            Mass Number = 2
            Relative Atomic Mass = 2.01410177812(12)
            Isotopic Composition = 0.000115(70)
            Standard Atomic Weight = [1.00784,1.00811]
            Notes = m

            Atomic Number = 1
            Atomic Symbol = T
            Mass Number = 3
            Relative Atomic Mass = 3.0160492779(24)
            Isotopic Composition =
            Standard Atomic Weight = [1.00784,1.00811]
            Notes = m

            Atomic Number = 2
            Atomic Symbol = He
            Mass Number = 3
            Relative Atomic Mass = 3.0160293201(25)
            Isotopic Composition = 0.00000134(3)
            Standard Atomic Weight = 4.002602(2)
            Notes = g,r

            Atomic Number = 2
            Atomic Symbol = He
            Mass Number = 4
            Relative Atomic Mass = 4.00260325413(6)
            Isotopic Composition = 0.99999866(3)
            Standard Atomic Weight = 4.002602(2)
            Notes = g,r

            Atomic Number = 3
            Atomic Symbol = Li
            Mass Number = 6
            Relative Atomic Mass = 6.0151228874(16)
            Isotopic Composition = 0.0759(4)
            Standard Atomic Weight = [6.938,6.997]
            Notes = m

            Atomic Number = 3
            Atomic Symbol = Li
            Mass Number = 7
            Relative Atomic Mass = 7.0160034366(45)
            Isotopic Composition = 0.9241(4)
            Standard Atomic Weight = [6.938,6.997]
            Notes = m

            Atomic Number = 4
            Atomic Symbol = Be
            Mass Number = 9
            Relative Atomic Mass = 9.012183065(82)
            Isotopic Composition = 1
            Standard Atomic Weight = 9.0121831(5)
            Notes = &nbsp;

            Atomic Number = 5
            Atomic Symbol = B
            Mass Number = 10
            Relative Atomic Mass = 10.01293695(41)
            Isotopic Composition = 0.199(7)
            Standard Atomic Weight = [10.806,10.821]
            Notes = m

            Atomic Number = 5
            Atomic Symbol = B
            Mass Number = 11
            Relative Atomic Mass = 11.00930536(45)
            Isotopic Composition = 0.801(7)
            Standard Atomic Weight = [10.806,10.821]
            Notes = m

            Atomic Number = 6
            Atomic Symbol = C
            Mass Number = 12
            Relative Atomic Mass = 12.0000000(00)
            Isotopic Composition = 0.9893(8)
            Standard Atomic Weight = [12.0096,12.0116]
            Notes = &nbsp;

            Atomic Number = 6
            Atomic Symbol = C
            Mass Number = 13
            Relative Atomic Mass = 13.00335483507(23)
            Isotopic Composition = 0.0107(8)
            Standard Atomic Weight = [12.0096,12.0116]
            Notes = &nbsp;

            Atomic Number = 6
            Atomic Symbol = C
            Mass Number = 14
            Relative Atomic Mass = 14.0032419884(40)
            Isotopic Composition =
            Standard Atomic Weight = [12.0096,12.0116]
            Notes = &nbsp;

            Atomic Number = 7
            Atomic Symbol = N
            Mass Number = 14
            Relative Atomic Mass = 14.00307400443(20)
            Isotopic Composition = 0.99636(20)
            Standard Atomic Weight = [14.00643,14.00728]
            Notes = &nbsp;

            Atomic Number = 7
            Atomic Symbol = N
            Mass Number = 15
            Relative Atomic Mass = 15.00010889888(64)
            Isotopic Composition = 0.00364(20)
            Standard Atomic Weight = [14.00643,14.00728]
            Notes = &nbsp;

            Atomic Number = 8
            Atomic Symbol = O
            Mass Number = 16
            Relative Atomic Mass = 15.99491461957(17)
            Isotopic Composition = 0.99757(16)
            Standard Atomic Weight = [15.99903,15.99977]
            Notes = &nbsp;

            Atomic Number = 8
            Atomic Symbol = O
            Mass Number = 17
            Relative Atomic Mass = 16.99913175650(69)
            Isotopic Composition = 0.00038(1)
            Standard Atomic Weight = [15.99903,15.99977]
            Notes = &nbsp;

            Atomic Number = 8
            Atomic Symbol = O
            Mass Number = 18
            Relative Atomic Mass = 17.99915961286(76)
            Isotopic Composition = 0.00205(14)
            Standard Atomic Weight = [15.99903,15.99977]
            Notes = &nbsp;

            Atomic Number = 9
            Atomic Symbol = F
            Mass Number = 19
            Relative Atomic Mass = 18.99840316273(92)
            Isotopic Composition = 1
            Standard Atomic Weight = 18.998403163(6)
            Notes = &nbsp;

            Atomic Number = 10
            Atomic Symbol = Ne
            Mass Number = 20
            Relative Atomic Mass = 19.9924401762(17)
            Isotopic Composition = 0.9048(3)
            Standard Atomic Weight = 20.1797(6)
            Notes = g,m

            Atomic Number = 10
            Atomic Symbol = Ne
            Mass Number = 21
            Relative Atomic Mass = 20.993846685(41)
            Isotopic Composition = 0.0027(1)
            Standard Atomic Weight = 20.1797(6)
            Notes = g,m

            Atomic Number = 10
            Atomic Symbol = Ne
            Mass Number = 22
            Relative Atomic Mass = 21.991385114(18)
            Isotopic Composition = 0.0925(3)
            Standard Atomic Weight = 20.1797(6)
            Notes = g,m

            Atomic Number = 11
            Atomic Symbol = Na
            Mass Number = 23
            Relative Atomic Mass = 22.9897692820(19)
            Isotopic Composition = 1
            Standard Atomic Weight = 22.98976928(2)
            Notes = &nbsp;

            Atomic Number = 12
            Atomic Symbol = Mg
            Mass Number = 24
            Relative Atomic Mass = 23.985041697(14)
            Isotopic Composition = 0.7899(4)
            Standard Atomic Weight = [24.304,24.307]
            Notes = &nbsp;

            Atomic Number = 12
            Atomic Symbol = Mg
            Mass Number = 25
            Relative Atomic Mass = 24.985836976(50)
            Isotopic Composition = 0.1000(1)
            Standard Atomic Weight = [24.304,24.307]
            Notes = &nbsp;

            Atomic Number = 12
            Atomic Symbol = Mg
            Mass Number = 26
            Relative Atomic Mass = 25.982592968(31)
            Isotopic Composition = 0.1101(3)
            Standard Atomic Weight = [24.304,24.307]
            Notes = &nbsp;

            Atomic Number = 13
            Atomic Symbol = Al
            Mass Number = 27
            Relative Atomic Mass = 26.98153853(11)
            Isotopic Composition = 1
            Standard Atomic Weight = 26.9815385(7)
            Notes = &nbsp;

            Atomic Number = 14
            Atomic Symbol = Si
            Mass Number = 28
            Relative Atomic Mass = 27.97692653465(44)
            Isotopic Composition = 0.92223(19)
            Standard Atomic Weight = [28.084,28.086]
            Notes = &nbsp;

            Atomic Number = 14
            Atomic Symbol = Si
            Mass Number = 29
            Relative Atomic Mass = 28.97649466490(52)
            Isotopic Composition = 0.04685(8)
            Standard Atomic Weight = [28.084,28.086]
            Notes = &nbsp;

            Atomic Number = 14
            Atomic Symbol = Si
            Mass Number = 30
            Relative Atomic Mass = 29.973770136(23)
            Isotopic Composition = 0.03092(11)
            Standard Atomic Weight = [28.084,28.086]
            Notes = &nbsp;

            Atomic Number = 15
            Atomic Symbol = P
            Mass Number = 31
            Relative Atomic Mass = 30.97376199842(70)
            Isotopic Composition = 1
            Standard Atomic Weight = 30.973761998(5)
            Notes = &nbsp;

            Atomic Number = 16
            Atomic Symbol = S
            Mass Number = 32
            Relative Atomic Mass = 31.9720711744(14)
            Isotopic Composition = 0.9499(26)
            Standard Atomic Weight = [32.059,32.076]
            Notes = &nbsp;

            Atomic Number = 16
            Atomic Symbol = S
            Mass Number = 33
            Relative Atomic Mass = 32.9714589098(15)
            Isotopic Composition = 0.0075(2)
            Standard Atomic Weight = [32.059,32.076]
            Notes = &nbsp;

            Atomic Number = 16
            Atomic Symbol = S
            Mass Number = 34
            Relative Atomic Mass = 33.967867004(47)
            Isotopic Composition = 0.0425(24)
            Standard Atomic Weight = [32.059,32.076]
            Notes = &nbsp;

            Atomic Number = 16
            Atomic Symbol = S
            Mass Number = 36
            Relative Atomic Mass = 35.96708071(20)
            Isotopic Composition = 0.0001(1)
            Standard Atomic Weight = [32.059,32.076]
            Notes = &nbsp;

            Atomic Number = 17
            Atomic Symbol = Cl
            Mass Number = 35
            Relative Atomic Mass = 34.968852682(37)
            Isotopic Composition = 0.7576(10)
            Standard Atomic Weight = [35.446,35.457]
            Notes = m

            Atomic Number = 17
            Atomic Symbol = Cl
            Mass Number = 37
            Relative Atomic Mass = 36.965902602(55)
            Isotopic Composition = 0.2424(10)
            Standard Atomic Weight = [35.446,35.457]
            Notes = m

            Atomic Number = 18
            Atomic Symbol = Ar
            Mass Number = 36
            Relative Atomic Mass = 35.967545105(28)
            Isotopic Composition = 0.003336(21)
            Standard Atomic Weight = 39.948(1)
            Notes = g,r

            Atomic Number = 18
            Atomic Symbol = Ar
            Mass Number = 38
            Relative Atomic Mass = 37.96273211(21)
            Isotopic Composition = 0.000629(7)
            Standard Atomic Weight = 39.948(1)
            Notes = g,r

            Atomic Number = 18
            Atomic Symbol = Ar
            Mass Number = 40
            Relative Atomic Mass = 39.9623831237(24)
            Isotopic Composition = 0.996035(25)
            Standard Atomic Weight = 39.948(1)
            Notes = g,r

            Atomic Number = 19
            Atomic Symbol = K
            Mass Number = 39
            Relative Atomic Mass = 38.9637064864(49)
            Isotopic Composition = 0.932581(44)
            Standard Atomic Weight = 39.0983(1)
            Notes = &nbsp;

            Atomic Number = 19
            Atomic Symbol = K
            Mass Number = 40
            Relative Atomic Mass = 39.963998166(60)
            Isotopic Composition = 0.000117(1)
            Standard Atomic Weight = 39.0983(1)
            Notes = &nbsp;

            Atomic Number = 19
            Atomic Symbol = K
            Mass Number = 41
            Relative Atomic Mass = 40.9618252579(41)
            Isotopic Composition = 0.067302(44)
            Standard Atomic Weight = 39.0983(1)
            Notes = &nbsp;

            Atomic Number = 20
            Atomic Symbol = Ca
            Mass Number = 40
            Relative Atomic Mass = 39.962590863(22)
            Isotopic Composition = 0.96941(156)
            Standard Atomic Weight = 40.078(4)
            Notes = g

            Atomic Number = 20
            Atomic Symbol = Ca
            Mass Number = 42
            Relative Atomic Mass = 41.95861783(16)
            Isotopic Composition = 0.00647(23)
            Standard Atomic Weight = 40.078(4)
            Notes = g

            Atomic Number = 20
            Atomic Symbol = Ca
            Mass Number = 43
            Relative Atomic Mass = 42.95876644(24)
            Isotopic Composition = 0.00135(10)
            Standard Atomic Weight = 40.078(4)
            Notes = g

            Atomic Number = 20
            Atomic Symbol = Ca
            Mass Number = 44
            Relative Atomic Mass = 43.95548156(35)
            Isotopic Composition = 0.02086(110)
            Standard Atomic Weight = 40.078(4)
            Notes = g

            Atomic Number = 20
            Atomic Symbol = Ca
            Mass Number = 46
            Relative Atomic Mass = 45.9536890(24)
            Isotopic Composition = 0.00004(3)
            Standard Atomic Weight = 40.078(4)
            Notes = g

            Atomic Number = 20
            Atomic Symbol = Ca
            Mass Number = 48
            Relative Atomic Mass = 47.95252276(13)
            Isotopic Composition = 0.00187(21)
            Standard Atomic Weight = 40.078(4)
            Notes = g

            Atomic Number = 21
            Atomic Symbol = Sc
            Mass Number = 45
            Relative Atomic Mass = 44.95590828(77)
            Isotopic Composition = 1
            Standard Atomic Weight = 44.955908(5)
            Notes = &nbsp;

            Atomic Number = 22
            Atomic Symbol = Ti
            Mass Number = 46
            Relative Atomic Mass = 45.95262772(35)
            Isotopic Composition = 0.0825(3)
            Standard Atomic Weight = 47.867(1)
            Notes = &nbsp;

            Atomic Number = 22
            Atomic Symbol = Ti
            Mass Number = 47
            Relative Atomic Mass = 46.95175879(38)
            Isotopic Composition = 0.0744(2)
            Standard Atomic Weight = 47.867(1)
            Notes = &nbsp;

            Atomic Number = 22
            Atomic Symbol = Ti
            Mass Number = 48
            Relative Atomic Mass = 47.94794198(38)
            Isotopic Composition = 0.7372(3)
            Standard Atomic Weight = 47.867(1)
            Notes = &nbsp;

            Atomic Number = 22
            Atomic Symbol = Ti
            Mass Number = 49
            Relative Atomic Mass = 48.94786568(39)
            Isotopic Composition = 0.0541(2)
            Standard Atomic Weight = 47.867(1)
            Notes = &nbsp;

            Atomic Number = 22
            Atomic Symbol = Ti
            Mass Number = 50
            Relative Atomic Mass = 49.94478689(39)
            Isotopic Composition = 0.0518(2)
            Standard Atomic Weight = 47.867(1)
            Notes = &nbsp;

            Atomic Number = 23
            Atomic Symbol = V
            Mass Number = 50
            Relative Atomic Mass = 49.94715601(95)
            Isotopic Composition = 0.00250(4)
            Standard Atomic Weight = 50.9415(1)
            Notes = &nbsp;

            Atomic Number = 23
            Atomic Symbol = V
            Mass Number = 51
            Relative Atomic Mass = 50.94395704(94)
            Isotopic Composition = 0.99750(4)
            Standard Atomic Weight = 50.9415(1)
            Notes = &nbsp;

            Atomic Number = 24
            Atomic Symbol = Cr
            Mass Number = 50
            Relative Atomic Mass = 49.94604183(94)
            Isotopic Composition = 0.04345(13)
            Standard Atomic Weight = 51.9961(6)
            Notes = &nbsp;

            Atomic Number = 24
            Atomic Symbol = Cr
            Mass Number = 52
            Relative Atomic Mass = 51.94050623(63)
            Isotopic Composition = 0.83789(18)
            Standard Atomic Weight = 51.9961(6)
            Notes = &nbsp;

            Atomic Number = 24
            Atomic Symbol = Cr
            Mass Number = 53
            Relative Atomic Mass = 52.94064815(62)
            Isotopic Composition = 0.09501(17)
            Standard Atomic Weight = 51.9961(6)
            Notes = &nbsp;

            Atomic Number = 24
            Atomic Symbol = Cr
            Mass Number = 54
            Relative Atomic Mass = 53.93887916(61)
            Isotopic Composition = 0.02365(7)
            Standard Atomic Weight = 51.9961(6)
            Notes = &nbsp;

            Atomic Number = 25
            Atomic Symbol = Mn
            Mass Number = 55
            Relative Atomic Mass = 54.93804391(48)
            Isotopic Composition = 1
            Standard Atomic Weight = 54.938044(3)
            Notes = &nbsp;

            Atomic Number = 26
            Atomic Symbol = Fe
            Mass Number = 54
            Relative Atomic Mass = 53.93960899(53)
            Isotopic Composition = 0.05845(35)
            Standard Atomic Weight = 55.845(2)
            Notes = &nbsp;

            Atomic Number = 26
            Atomic Symbol = Fe
            Mass Number = 56
            Relative Atomic Mass = 55.93493633(49)
            Isotopic Composition = 0.91754(36)
            Standard Atomic Weight = 55.845(2)
            Notes = &nbsp;

            Atomic Number = 26
            Atomic Symbol = Fe
            Mass Number = 57
            Relative Atomic Mass = 56.93539284(49)
            Isotopic Composition = 0.02119(10)
            Standard Atomic Weight = 55.845(2)
            Notes = &nbsp;

            Atomic Number = 26
            Atomic Symbol = Fe
            Mass Number = 58
            Relative Atomic Mass = 57.93327443(53)
            Isotopic Composition = 0.00282(4)
            Standard Atomic Weight = 55.845(2)
            Notes = &nbsp;

            Atomic Number = 27
            Atomic Symbol = Co
            Mass Number = 59
            Relative Atomic Mass = 58.93319429(56)
            Isotopic Composition = 1
            Standard Atomic Weight = 58.933194(4)
            Notes = &nbsp;

            Atomic Number = 28
            Atomic Symbol = Ni
            Mass Number = 58
            Relative Atomic Mass = 57.93534241(52)
            Isotopic Composition = 0.68077(19)
            Standard Atomic Weight = 58.6934(4)
            Notes = r

            Atomic Number = 28
            Atomic Symbol = Ni
            Mass Number = 60
            Relative Atomic Mass = 59.93078588(52)
            Isotopic Composition = 0.26223(15)
            Standard Atomic Weight = 58.6934(4)
            Notes = r

            Atomic Number = 28
            Atomic Symbol = Ni
            Mass Number = 61
            Relative Atomic Mass = 60.93105557(52)
            Isotopic Composition = 0.011399(13)
            Standard Atomic Weight = 58.6934(4)
            Notes = r

            Atomic Number = 28
            Atomic Symbol = Ni
            Mass Number = 62
            Relative Atomic Mass = 61.92834537(55)
            Isotopic Composition = 0.036346(40)
            Standard Atomic Weight = 58.6934(4)
            Notes = r

            Atomic Number = 28
            Atomic Symbol = Ni
            Mass Number = 64
            Relative Atomic Mass = 63.92796682(58)
            Isotopic Composition = 0.009255(19)
            Standard Atomic Weight = 58.6934(4)
            Notes = r

            Atomic Number = 29
            Atomic Symbol = Cu
            Mass Number = 63
            Relative Atomic Mass = 62.92959772(56)
            Isotopic Composition = 0.6915(15)
            Standard Atomic Weight = 63.546(3)
            Notes = r

            Atomic Number = 29
            Atomic Symbol = Cu
            Mass Number = 65
            Relative Atomic Mass = 64.92778970(71)
            Isotopic Composition = 0.3085(15)
            Standard Atomic Weight = 63.546(3)
            Notes = r

            Atomic Number = 30
            Atomic Symbol = Zn
            Mass Number = 64
            Relative Atomic Mass = 63.92914201(71)
            Isotopic Composition = 0.4917(75)
            Standard Atomic Weight = 65.38(2)
            Notes = r

            Atomic Number = 30
            Atomic Symbol = Zn
            Mass Number = 66
            Relative Atomic Mass = 65.92603381(94)
            Isotopic Composition = 0.2773(98)
            Standard Atomic Weight = 65.38(2)
            Notes = r

            Atomic Number = 30
            Atomic Symbol = Zn
            Mass Number = 67
            Relative Atomic Mass = 66.92712775(96)
            Isotopic Composition = 0.0404(16)
            Standard Atomic Weight = 65.38(2)
            Notes = r

            Atomic Number = 30
            Atomic Symbol = Zn
            Mass Number = 68
            Relative Atomic Mass = 67.92484455(98)
            Isotopic Composition = 0.1845(63)
            Standard Atomic Weight = 65.38(2)
            Notes = r

            Atomic Number = 30
            Atomic Symbol = Zn
            Mass Number = 70
            Relative Atomic Mass = 69.9253192(21)
            Isotopic Composition = 0.0061(10)
            Standard Atomic Weight = 65.38(2)
            Notes = r

            Atomic Number = 31
            Atomic Symbol = Ga
            Mass Number = 69
            Relative Atomic Mass = 68.9255735(13)
            Isotopic Composition = 0.60108(9)
            Standard Atomic Weight = 69.723(1)
            Notes = &nbsp;

            Atomic Number = 31
            Atomic Symbol = Ga
            Mass Number = 71
            Relative Atomic Mass = 70.92470258(87)
            Isotopic Composition = 0.39892(9)
            Standard Atomic Weight = 69.723(1)
            Notes = &nbsp;

            Atomic Number = 32
            Atomic Symbol = Ge
            Mass Number = 70
            Relative Atomic Mass = 69.92424875(90)
            Isotopic Composition = 0.2057(27)
            Standard Atomic Weight = 72.630(8)
            Notes = &nbsp;

            Atomic Number = 32
            Atomic Symbol = Ge
            Mass Number = 72
            Relative Atomic Mass = 71.922075826(81)
            Isotopic Composition = 0.2745(32)
            Standard Atomic Weight = 72.630(8)
            Notes = &nbsp;

            Atomic Number = 32
            Atomic Symbol = Ge
            Mass Number = 73
            Relative Atomic Mass = 72.923458956(61)
            Isotopic Composition = 0.0775(12)
            Standard Atomic Weight = 72.630(8)
            Notes = &nbsp;

            Atomic Number = 32
            Atomic Symbol = Ge
            Mass Number = 74
            Relative Atomic Mass = 73.921177761(13)
            Isotopic Composition = 0.3650(20)
            Standard Atomic Weight = 72.630(8)
            Notes = &nbsp;

            Atomic Number = 32
            Atomic Symbol = Ge
            Mass Number = 76
            Relative Atomic Mass = 75.921402726(19)
            Isotopic Composition = 0.0773(12)
            Standard Atomic Weight = 72.630(8)
            Notes = &nbsp;

            Atomic Number = 33
            Atomic Symbol = As
            Mass Number = 75
            Relative Atomic Mass = 74.92159457(95)
            Isotopic Composition = 1
            Standard Atomic Weight = 74.921595(6)
            Notes = &nbsp;

            Atomic Number = 34
            Atomic Symbol = Se
            Mass Number = 74
            Relative Atomic Mass = 73.922475934(15)
            Isotopic Composition = 0.0089(4)
            Standard Atomic Weight = 78.971(8)
            Notes = r

            Atomic Number = 34
            Atomic Symbol = Se
            Mass Number = 76
            Relative Atomic Mass = 75.919213704(17)
            Isotopic Composition = 0.0937(29)
            Standard Atomic Weight = 78.971(8)
            Notes = r

            Atomic Number = 34
            Atomic Symbol = Se
            Mass Number = 77
            Relative Atomic Mass = 76.919914154(67)
            Isotopic Composition = 0.0763(16)
            Standard Atomic Weight = 78.971(8)
            Notes = r

            Atomic Number = 34
            Atomic Symbol = Se
            Mass Number = 78
            Relative Atomic Mass = 77.91730928(20)
            Isotopic Composition = 0.2377(28)
            Standard Atomic Weight = 78.971(8)
            Notes = r

            Atomic Number = 34
            Atomic Symbol = Se
            Mass Number = 80
            Relative Atomic Mass = 79.9165218(13)
            Isotopic Composition = 0.4961(41)
            Standard Atomic Weight = 78.971(8)
            Notes = r

            Atomic Number = 34
            Atomic Symbol = Se
            Mass Number = 82
            Relative Atomic Mass = 81.9166995(15)
            Isotopic Composition = 0.0873(22)
            Standard Atomic Weight = 78.971(8)
            Notes = r

            Atomic Number = 35
            Atomic Symbol = Br
            Mass Number = 79
            Relative Atomic Mass = 78.9183376(14)
            Isotopic Composition = 0.5069(7)
            Standard Atomic Weight = [79.901,79.907]
            Notes = &nbsp;

            Atomic Number = 35
            Atomic Symbol = Br
            Mass Number = 81
            Relative Atomic Mass = 80.9162897(14)
            Isotopic Composition = 0.4931(7)
            Standard Atomic Weight = [79.901,79.907]
            Notes = &nbsp;

            Atomic Number = 36
            Atomic Symbol = Kr
            Mass Number = 78
            Relative Atomic Mass = 77.92036494(76)
            Isotopic Composition = 0.00355(3)
            Standard Atomic Weight = 83.798(2)
            Notes = g,m

            Atomic Number = 36
            Atomic Symbol = Kr
            Mass Number = 80
            Relative Atomic Mass = 79.91637808(75)
            Isotopic Composition = 0.02286(10)
            Standard Atomic Weight = 83.798(2)
            Notes = g,m

            Atomic Number = 36
            Atomic Symbol = Kr
            Mass Number = 82
            Relative Atomic Mass = 81.91348273(94)
            Isotopic Composition = 0.11593(31)
            Standard Atomic Weight = 83.798(2)
            Notes = g,m

            Atomic Number = 36
            Atomic Symbol = Kr
            Mass Number = 83
            Relative Atomic Mass = 82.91412716(32)
            Isotopic Composition = 0.11500(19)
            Standard Atomic Weight = 83.798(2)
            Notes = g,m

            Atomic Number = 36
            Atomic Symbol = Kr
            Mass Number = 84
            Relative Atomic Mass = 83.9114977282(44)
            Isotopic Composition = 0.56987(15)
            Standard Atomic Weight = 83.798(2)
            Notes = g,m

            Atomic Number = 36
            Atomic Symbol = Kr
            Mass Number = 86
            Relative Atomic Mass = 85.9106106269(41)
            Isotopic Composition = 0.17279(41)
            Standard Atomic Weight = 83.798(2)
            Notes = g,m

            Atomic Number = 37
            Atomic Symbol = Rb
            Mass Number = 85
            Relative Atomic Mass = 84.9117897379(54)
            Isotopic Composition = 0.7217(2)
            Standard Atomic Weight = 85.4678(3)
            Notes = g

            Atomic Number = 37
            Atomic Symbol = Rb
            Mass Number = 87
            Relative Atomic Mass = 86.9091805310(60)
            Isotopic Composition = 0.2783(2)
            Standard Atomic Weight = 85.4678(3)
            Notes = g

            Atomic Number = 38
            Atomic Symbol = Sr
            Mass Number = 84
            Relative Atomic Mass = 83.9134191(13)
            Isotopic Composition = 0.0056(1)
            Standard Atomic Weight = 87.62(1)
            Notes = g,r

            Atomic Number = 38
            Atomic Symbol = Sr
            Mass Number = 86
            Relative Atomic Mass = 85.9092606(12)
            Isotopic Composition = 0.0986(1)
            Standard Atomic Weight = 87.62(1)
            Notes = g,r

            Atomic Number = 38
            Atomic Symbol = Sr
            Mass Number = 87
            Relative Atomic Mass = 86.9088775(12)
            Isotopic Composition = 0.0700(1)
            Standard Atomic Weight = 87.62(1)
            Notes = g,r

            Atomic Number = 38
            Atomic Symbol = Sr
            Mass Number = 88
            Relative Atomic Mass = 87.9056125(12)
            Isotopic Composition = 0.8258(1)
            Standard Atomic Weight = 87.62(1)
            Notes = g,r

            Atomic Number = 39
            Atomic Symbol = Y
            Mass Number = 89
            Relative Atomic Mass = 88.9058403(24)
            Isotopic Composition = 1
            Standard Atomic Weight = 88.90584(2)
            Notes = &nbsp;

            Atomic Number = 40
            Atomic Symbol = Zr
            Mass Number = 90
            Relative Atomic Mass = 89.9046977(20)
            Isotopic Composition = 0.5145(40)
            Standard Atomic Weight = 91.224(2)
            Notes = g

            Atomic Number = 40
            Atomic Symbol = Zr
            Mass Number = 91
            Relative Atomic Mass = 90.9056396(20)
            Isotopic Composition = 0.1122(5)
            Standard Atomic Weight = 91.224(2)
            Notes = g

            Atomic Number = 40
            Atomic Symbol = Zr
            Mass Number = 92
            Relative Atomic Mass = 91.9050347(20)
            Isotopic Composition = 0.1715(8)
            Standard Atomic Weight = 91.224(2)
            Notes = g

            Atomic Number = 40
            Atomic Symbol = Zr
            Mass Number = 94
            Relative Atomic Mass = 93.9063108(20)
            Isotopic Composition = 0.1738(28)
            Standard Atomic Weight = 91.224(2)
            Notes = g

            Atomic Number = 40
            Atomic Symbol = Zr
            Mass Number = 96
            Relative Atomic Mass = 95.9082714(21)
            Isotopic Composition = 0.0280(9)
            Standard Atomic Weight = 91.224(2)
            Notes = g

            Atomic Number = 41
            Atomic Symbol = Nb
            Mass Number = 93
            Relative Atomic Mass = 92.9063730(20)
            Isotopic Composition = 1
            Standard Atomic Weight = 92.90637(2)
            Notes = &nbsp;

            Atomic Number = 42
            Atomic Symbol = Mo
            Mass Number = 92
            Relative Atomic Mass = 91.90680796(84)
            Isotopic Composition = 0.1453(30)
            Standard Atomic Weight = 95.95(1)
            Notes = g

            Atomic Number = 42
            Atomic Symbol = Mo
            Mass Number = 94
            Relative Atomic Mass = 93.90508490(48)
            Isotopic Composition = 0.0915(9)
            Standard Atomic Weight = 95.95(1)
            Notes = g

            Atomic Number = 42
            Atomic Symbol = Mo
            Mass Number = 95
            Relative Atomic Mass = 94.90583877(47)
            Isotopic Composition = 0.1584(11)
            Standard Atomic Weight = 95.95(1)
            Notes = g

            Atomic Number = 42
            Atomic Symbol = Mo
            Mass Number = 96
            Relative Atomic Mass = 95.90467612(47)
            Isotopic Composition = 0.1667(15)
            Standard Atomic Weight = 95.95(1)
            Notes = g

            Atomic Number = 42
            Atomic Symbol = Mo
            Mass Number = 97
            Relative Atomic Mass = 96.90601812(49)
            Isotopic Composition = 0.0960(14)
            Standard Atomic Weight = 95.95(1)
            Notes = g

            Atomic Number = 42
            Atomic Symbol = Mo
            Mass Number = 98
            Relative Atomic Mass = 97.90540482(49)
            Isotopic Composition = 0.2439(37)
            Standard Atomic Weight = 95.95(1)
            Notes = g

            Atomic Number = 42
            Atomic Symbol = Mo
            Mass Number = 100
            Relative Atomic Mass = 99.9074718(11)
            Isotopic Composition = 0.0982(31)
            Standard Atomic Weight = 95.95(1)
            Notes = g

            Atomic Number = 43
            Atomic Symbol = Tc
            Mass Number = 97
            Relative Atomic Mass = 96.9063667(40)
            Isotopic Composition =
            Standard Atomic Weight = [98]
            Notes = &nbsp;

            Atomic Number = 43
            Atomic Symbol = Tc
            Mass Number = 98
            Relative Atomic Mass = 97.9072124(36)
            Isotopic Composition =
            Standard Atomic Weight = [98]
            Notes = &nbsp;

            Atomic Number = 43
            Atomic Symbol = Tc
            Mass Number = 99
            Relative Atomic Mass = 98.9062508(10)
            Isotopic Composition =
            Standard Atomic Weight = [98]
            Notes = &nbsp;

            Atomic Number = 44
            Atomic Symbol = Ru
            Mass Number = 96
            Relative Atomic Mass = 95.90759025(49)
            Isotopic Composition = 0.0554(14)
            Standard Atomic Weight = 101.07(2)
            Notes = g

            Atomic Number = 44
            Atomic Symbol = Ru
            Mass Number = 98
            Relative Atomic Mass = 97.9052868(69)
            Isotopic Composition = 0.0187(3)
            Standard Atomic Weight = 101.07(2)
            Notes = g

            Atomic Number = 44
            Atomic Symbol = Ru
            Mass Number = 99
            Relative Atomic Mass = 98.9059341(11)
            Isotopic Composition = 0.1276(14)
            Standard Atomic Weight = 101.07(2)
            Notes = g

            Atomic Number = 44
            Atomic Symbol = Ru
            Mass Number = 100
            Relative Atomic Mass = 99.9042143(11)
            Isotopic Composition = 0.1260(7)
            Standard Atomic Weight = 101.07(2)
            Notes = g

            Atomic Number = 44
            Atomic Symbol = Ru
            Mass Number = 101
            Relative Atomic Mass = 100.9055769(12)
            Isotopic Composition = 0.1706(2)
            Standard Atomic Weight = 101.07(2)
            Notes = g

            Atomic Number = 44
            Atomic Symbol = Ru
            Mass Number = 102
            Relative Atomic Mass = 101.9043441(12)
            Isotopic Composition = 0.3155(14)
            Standard Atomic Weight = 101.07(2)
            Notes = g

            Atomic Number = 44
            Atomic Symbol = Ru
            Mass Number = 104
            Relative Atomic Mass = 103.9054275(28)
            Isotopic Composition = 0.1862(27)
            Standard Atomic Weight = 101.07(2)
            Notes = g

            Atomic Number = 45
            Atomic Symbol = Rh
            Mass Number = 103
            Relative Atomic Mass = 102.9054980(26)
            Isotopic Composition = 1
            Standard Atomic Weight = 102.90550(2)
            Notes = &nbsp;

            Atomic Number = 46
            Atomic Symbol = Pd
            Mass Number = 102
            Relative Atomic Mass = 101.9056022(28)
            Isotopic Composition = 0.0102(1)
            Standard Atomic Weight = 106.42(1)
            Notes = g

            Atomic Number = 46
            Atomic Symbol = Pd
            Mass Number = 104
            Relative Atomic Mass = 103.9040305(14)
            Isotopic Composition = 0.1114(8)
            Standard Atomic Weight = 106.42(1)
            Notes = g

            Atomic Number = 46
            Atomic Symbol = Pd
            Mass Number = 105
            Relative Atomic Mass = 104.9050796(12)
            Isotopic Composition = 0.2233(8)
            Standard Atomic Weight = 106.42(1)
            Notes = g

            Atomic Number = 46
            Atomic Symbol = Pd
            Mass Number = 106
            Relative Atomic Mass = 105.9034804(12)
            Isotopic Composition = 0.2733(3)
            Standard Atomic Weight = 106.42(1)
            Notes = g

            Atomic Number = 46
            Atomic Symbol = Pd
            Mass Number = 108
            Relative Atomic Mass = 107.9038916(12)
            Isotopic Composition = 0.2646(9)
            Standard Atomic Weight = 106.42(1)
            Notes = g

            Atomic Number = 46
            Atomic Symbol = Pd
            Mass Number = 110
            Relative Atomic Mass = 109.90517220(75)
            Isotopic Composition = 0.1172(9)
            Standard Atomic Weight = 106.42(1)
            Notes = g

            Atomic Number = 47
            Atomic Symbol = Ag
            Mass Number = 107
            Relative Atomic Mass = 106.9050916(26)
            Isotopic Composition = 0.51839(8)
            Standard Atomic Weight = 107.8682(2)
            Notes = g

            Atomic Number = 47
            Atomic Symbol = Ag
            Mass Number = 109
            Relative Atomic Mass = 108.9047553(14)
            Isotopic Composition = 0.48161(8)
            Standard Atomic Weight = 107.8682(2)
            Notes = g

            Atomic Number = 48
            Atomic Symbol = Cd
            Mass Number = 106
            Relative Atomic Mass = 105.9064599(12)
            Isotopic Composition = 0.0125(6)
            Standard Atomic Weight = 112.414(4)
            Notes = g

            Atomic Number = 48
            Atomic Symbol = Cd
            Mass Number = 108
            Relative Atomic Mass = 107.9041834(12)
            Isotopic Composition = 0.0089(3)
            Standard Atomic Weight = 112.414(4)
            Notes = g

            Atomic Number = 48
            Atomic Symbol = Cd
            Mass Number = 110
            Relative Atomic Mass = 109.90300661(61)
            Isotopic Composition = 0.1249(18)
            Standard Atomic Weight = 112.414(4)
            Notes = g

            Atomic Number = 48
            Atomic Symbol = Cd
            Mass Number = 111
            Relative Atomic Mass = 110.90418287(61)
            Isotopic Composition = 0.1280(12)
            Standard Atomic Weight = 112.414(4)
            Notes = g

            Atomic Number = 48
            Atomic Symbol = Cd
            Mass Number = 112
            Relative Atomic Mass = 111.90276287(60)
            Isotopic Composition = 0.2413(21)
            Standard Atomic Weight = 112.414(4)
            Notes = g

            Atomic Number = 48
            Atomic Symbol = Cd
            Mass Number = 113
            Relative Atomic Mass = 112.90440813(45)
            Isotopic Composition = 0.1222(12)
            Standard Atomic Weight = 112.414(4)
            Notes = g

            Atomic Number = 48
            Atomic Symbol = Cd
            Mass Number = 114
            Relative Atomic Mass = 113.90336509(43)
            Isotopic Composition = 0.2873(42)
            Standard Atomic Weight = 112.414(4)
            Notes = g

            Atomic Number = 48
            Atomic Symbol = Cd
            Mass Number = 116
            Relative Atomic Mass = 115.90476315(17)
            Isotopic Composition = 0.0749(18)
            Standard Atomic Weight = 112.414(4)
            Notes = g

            Atomic Number = 49
            Atomic Symbol = In
            Mass Number = 113
            Relative Atomic Mass = 112.90406184(91)
            Isotopic Composition = 0.0429(5)
            Standard Atomic Weight = 114.818(1)
            Notes = &nbsp;

            Atomic Number = 49
            Atomic Symbol = In
            Mass Number = 115
            Relative Atomic Mass = 114.903878776(12)
            Isotopic Composition = 0.9571(5)
            Standard Atomic Weight = 114.818(1)
            Notes = &nbsp;

            Atomic Number = 50
            Atomic Symbol = Sn
            Mass Number = 112
            Relative Atomic Mass = 111.90482387(61)
            Isotopic Composition = 0.0097(1)
            Standard Atomic Weight = 118.710(7)
            Notes = g

            Atomic Number = 50
            Atomic Symbol = Sn
            Mass Number = 114
            Relative Atomic Mass = 113.9027827(10)
            Isotopic Composition = 0.0066(1)
            Standard Atomic Weight = 118.710(7)
            Notes = g

            Atomic Number = 50
            Atomic Symbol = Sn
            Mass Number = 115
            Relative Atomic Mass = 114.903344699(16)
            Isotopic Composition = 0.0034(1)
            Standard Atomic Weight = 118.710(7)
            Notes = g

            Atomic Number = 50
            Atomic Symbol = Sn
            Mass Number = 116
            Relative Atomic Mass = 115.90174280(10)
            Isotopic Composition = 0.1454(9)
            Standard Atomic Weight = 118.710(7)
            Notes = g

            Atomic Number = 50
            Atomic Symbol = Sn
            Mass Number = 117
            Relative Atomic Mass = 116.90295398(52)
            Isotopic Composition = 0.0768(7)
            Standard Atomic Weight = 118.710(7)
            Notes = g

            Atomic Number = 50
            Atomic Symbol = Sn
            Mass Number = 118
            Relative Atomic Mass = 117.90160657(54)
            Isotopic Composition = 0.2422(9)
            Standard Atomic Weight = 118.710(7)
            Notes = g

            Atomic Number = 50
            Atomic Symbol = Sn
            Mass Number = 119
            Relative Atomic Mass = 118.90331117(78)
            Isotopic Composition = 0.0859(4)
            Standard Atomic Weight = 118.710(7)
            Notes = g

            Atomic Number = 50
            Atomic Symbol = Sn
            Mass Number = 120
            Relative Atomic Mass = 119.90220163(97)
            Isotopic Composition = 0.3258(9)
            Standard Atomic Weight = 118.710(7)
            Notes = g

            Atomic Number = 50
            Atomic Symbol = Sn
            Mass Number = 122
            Relative Atomic Mass = 121.9034438(26)
            Isotopic Composition = 0.0463(3)
            Standard Atomic Weight = 118.710(7)
            Notes = g

            Atomic Number = 50
            Atomic Symbol = Sn
            Mass Number = 124
            Relative Atomic Mass = 123.9052766(11)
            Isotopic Composition = 0.0579(5)
            Standard Atomic Weight = 118.710(7)
            Notes = g

            Atomic Number = 51
            Atomic Symbol = Sb
            Mass Number = 121
            Relative Atomic Mass = 120.9038120(30)
            Isotopic Composition = 0.5721(5)
            Standard Atomic Weight = 121.760(1)
            Notes = g

            Atomic Number = 51
            Atomic Symbol = Sb
            Mass Number = 123
            Relative Atomic Mass = 122.9042132(23)
            Isotopic Composition = 0.4279(5)
            Standard Atomic Weight = 121.760(1)
            Notes = g

            Atomic Number = 52
            Atomic Symbol = Te
            Mass Number = 120
            Relative Atomic Mass = 119.9040593(33)
            Isotopic Composition = 0.0009(1)
            Standard Atomic Weight = 127.60(3)
            Notes = g

            Atomic Number = 52
            Atomic Symbol = Te
            Mass Number = 122
            Relative Atomic Mass = 121.9030435(16)
            Isotopic Composition = 0.0255(12)
            Standard Atomic Weight = 127.60(3)
            Notes = g

            Atomic Number = 52
            Atomic Symbol = Te
            Mass Number = 123
            Relative Atomic Mass = 122.9042698(16)
            Isotopic Composition = 0.0089(3)
            Standard Atomic Weight = 127.60(3)
            Notes = g

            Atomic Number = 52
            Atomic Symbol = Te
            Mass Number = 124
            Relative Atomic Mass = 123.9028171(16)
            Isotopic Composition = 0.0474(14)
            Standard Atomic Weight = 127.60(3)
            Notes = g

            Atomic Number = 52
            Atomic Symbol = Te
            Mass Number = 125
            Relative Atomic Mass = 124.9044299(16)
            Isotopic Composition = 0.0707(15)
            Standard Atomic Weight = 127.60(3)
            Notes = g

            Atomic Number = 52
            Atomic Symbol = Te
            Mass Number = 126
            Relative Atomic Mass = 125.9033109(16)
            Isotopic Composition = 0.1884(25)
            Standard Atomic Weight = 127.60(3)
            Notes = g

            Atomic Number = 52
            Atomic Symbol = Te
            Mass Number = 128
            Relative Atomic Mass = 127.90446128(93)
            Isotopic Composition = 0.3174(8)
            Standard Atomic Weight = 127.60(3)
            Notes = g

            Atomic Number = 52
            Atomic Symbol = Te
            Mass Number = 130
            Relative Atomic Mass = 129.906222748(12)
            Isotopic Composition = 0.3408(62)
            Standard Atomic Weight = 127.60(3)
            Notes = g

            Atomic Number = 53
            Atomic Symbol = I
            Mass Number = 127
            Relative Atomic Mass = 126.9044719(39)
            Isotopic Composition = 1
            Standard Atomic Weight = 126.90447(3)
            Notes = &nbsp;

            Atomic Number = 54
            Atomic Symbol = Xe
            Mass Number = 124
            Relative Atomic Mass = 123.9058920(19)
            Isotopic Composition = 0.000952(3)
            Standard Atomic Weight = 131.293(6)
            Notes = g,m

            Atomic Number = 54
            Atomic Symbol = Xe
            Mass Number = 126
            Relative Atomic Mass = 125.9042983(38)
            Isotopic Composition = 0.000890(2)
            Standard Atomic Weight = 131.293(6)
            Notes = g,m

            Atomic Number = 54
            Atomic Symbol = Xe
            Mass Number = 128
            Relative Atomic Mass = 127.9035310(11)
            Isotopic Composition = 0.019102(8)
            Standard Atomic Weight = 131.293(6)
            Notes = g,m

            Atomic Number = 54
            Atomic Symbol = Xe
            Mass Number = 129
            Relative Atomic Mass = 128.9047808611(60)
            Isotopic Composition = 0.264006(82)
            Standard Atomic Weight = 131.293(6)
            Notes = g,m

            Atomic Number = 54
            Atomic Symbol = Xe
            Mass Number = 130
            Relative Atomic Mass = 129.903509349(10)
            Isotopic Composition = 0.040710(13)
            Standard Atomic Weight = 131.293(6)
            Notes = g,m

            Atomic Number = 54
            Atomic Symbol = Xe
            Mass Number = 131
            Relative Atomic Mass = 130.90508406(24)
            Isotopic Composition = 0.212324(30)
            Standard Atomic Weight = 131.293(6)
            Notes = g,m

            Atomic Number = 54
            Atomic Symbol = Xe
            Mass Number = 132
            Relative Atomic Mass = 131.9041550856(56)
            Isotopic Composition = 0.269086(33)
            Standard Atomic Weight = 131.293(6)
            Notes = g,m

            Atomic Number = 54
            Atomic Symbol = Xe
            Mass Number = 134
            Relative Atomic Mass = 133.90539466(90)
            Isotopic Composition = 0.104357(21)
            Standard Atomic Weight = 131.293(6)
            Notes = g,m

            Atomic Number = 54
            Atomic Symbol = Xe
            Mass Number = 136
            Relative Atomic Mass = 135.907214484(11)
            Isotopic Composition = 0.088573(44)
            Standard Atomic Weight = 131.293(6)
            Notes = g,m

            Atomic Number = 55
            Atomic Symbol = Cs
            Mass Number = 133
            Relative Atomic Mass = 132.9054519610(80)
            Isotopic Composition = 1
            Standard Atomic Weight = 132.90545196(6)
            Notes = &nbsp;

            Atomic Number = 56
            Atomic Symbol = Ba
            Mass Number = 130
            Relative Atomic Mass = 129.9063207(28)
            Isotopic Composition = 0.00106(1)
            Standard Atomic Weight = 137.327(7)
            Notes = &nbsp;

            Atomic Number = 56
            Atomic Symbol = Ba
            Mass Number = 132
            Relative Atomic Mass = 131.9050611(11)
            Isotopic Composition = 0.00101(1)
            Standard Atomic Weight = 137.327(7)
            Notes = &nbsp;

            Atomic Number = 56
            Atomic Symbol = Ba
            Mass Number = 134
            Relative Atomic Mass = 133.90450818(30)
            Isotopic Composition = 0.02417(18)
            Standard Atomic Weight = 137.327(7)
            Notes = &nbsp;

            Atomic Number = 56
            Atomic Symbol = Ba
            Mass Number = 135
            Relative Atomic Mass = 134.90568838(29)
            Isotopic Composition = 0.06592(12)
            Standard Atomic Weight = 137.327(7)
            Notes = &nbsp;

            Atomic Number = 56
            Atomic Symbol = Ba
            Mass Number = 136
            Relative Atomic Mass = 135.90457573(29)
            Isotopic Composition = 0.07854(24)
            Standard Atomic Weight = 137.327(7)
            Notes = &nbsp;

            Atomic Number = 56
            Atomic Symbol = Ba
            Mass Number = 137
            Relative Atomic Mass = 136.90582714(30)
            Isotopic Composition = 0.11232(24)
            Standard Atomic Weight = 137.327(7)
            Notes = &nbsp;

            Atomic Number = 56
            Atomic Symbol = Ba
            Mass Number = 138
            Relative Atomic Mass = 137.90524700(31)
            Isotopic Composition = 0.71698(42)
            Standard Atomic Weight = 137.327(7)
            Notes = &nbsp;

            Atomic Number = 57
            Atomic Symbol = La
            Mass Number = 138
            Relative Atomic Mass = 137.9071149(37)
            Isotopic Composition = 0.0008881(71)
            Standard Atomic Weight = 138.90547(7)
            Notes = g

            Atomic Number = 57
            Atomic Symbol = La
            Mass Number = 139
            Relative Atomic Mass = 138.9063563(24)
            Isotopic Composition = 0.9991119(71)
            Standard Atomic Weight = 138.90547(7)
            Notes = g

            Atomic Number = 58
            Atomic Symbol = Ce
            Mass Number = 136
            Relative Atomic Mass = 135.90712921(41)
            Isotopic Composition = 0.00185(2)
            Standard Atomic Weight = 140.116(1)
            Notes = g

            Atomic Number = 58
            Atomic Symbol = Ce
            Mass Number = 138
            Relative Atomic Mass = 137.905991(11)
            Isotopic Composition = 0.00251(2)
            Standard Atomic Weight = 140.116(1)
            Notes = g

            Atomic Number = 58
            Atomic Symbol = Ce
            Mass Number = 140
            Relative Atomic Mass = 139.9054431(23)
            Isotopic Composition = 0.88450(51)
            Standard Atomic Weight = 140.116(1)
            Notes = g

            Atomic Number = 58
            Atomic Symbol = Ce
            Mass Number = 142
            Relative Atomic Mass = 141.9092504(29)
            Isotopic Composition = 0.11114(51)
            Standard Atomic Weight = 140.116(1)
            Notes = g

            Atomic Number = 59
            Atomic Symbol = Pr
            Mass Number = 141
            Relative Atomic Mass = 140.9076576(23)
            Isotopic Composition = 1
            Standard Atomic Weight = 140.90766(2)
            Notes = &nbsp;

            Atomic Number = 60
            Atomic Symbol = Nd
            Mass Number = 142
            Relative Atomic Mass = 141.9077290(20)
            Isotopic Composition = 0.27152(40)
            Standard Atomic Weight = 144.242(3)
            Notes = g

            Atomic Number = 60
            Atomic Symbol = Nd
            Mass Number = 143
            Relative Atomic Mass = 142.9098200(20)
            Isotopic Composition = 0.12174(26)
            Standard Atomic Weight = 144.242(3)
            Notes = g

            Atomic Number = 60
            Atomic Symbol = Nd
            Mass Number = 144
            Relative Atomic Mass = 143.9100930(20)
            Isotopic Composition = 0.23798(19)
            Standard Atomic Weight = 144.242(3)
            Notes = g

            Atomic Number = 60
            Atomic Symbol = Nd
            Mass Number = 145
            Relative Atomic Mass = 144.9125793(20)
            Isotopic Composition = 0.08293(12)
            Standard Atomic Weight = 144.242(3)
            Notes = g

            Atomic Number = 60
            Atomic Symbol = Nd
            Mass Number = 146
            Relative Atomic Mass = 145.9131226(20)
            Isotopic Composition = 0.17189(32)
            Standard Atomic Weight = 144.242(3)
            Notes = g

            Atomic Number = 60
            Atomic Symbol = Nd
            Mass Number = 148
            Relative Atomic Mass = 147.9168993(26)
            Isotopic Composition = 0.05756(21)
            Standard Atomic Weight = 144.242(3)
            Notes = g

            Atomic Number = 60
            Atomic Symbol = Nd
            Mass Number = 150
            Relative Atomic Mass = 149.9209022(18)
            Isotopic Composition = 0.05638(28)
            Standard Atomic Weight = 144.242(3)
            Notes = g

            Atomic Number = 61
            Atomic Symbol = Pm
            Mass Number = 145
            Relative Atomic Mass = 144.9127559(33)
            Isotopic Composition =
            Standard Atomic Weight = [145]
            Notes = &nbsp;

            Atomic Number = 61
            Atomic Symbol = Pm
            Mass Number = 147
            Relative Atomic Mass = 146.9151450(19)
            Isotopic Composition =
            Standard Atomic Weight = [145]
            Notes = &nbsp;

            Atomic Number = 62
            Atomic Symbol = Sm
            Mass Number = 144
            Relative Atomic Mass = 143.9120065(21)
            Isotopic Composition = 0.0307(7)
            Standard Atomic Weight = 150.36(2)
            Notes = g

            Atomic Number = 62
            Atomic Symbol = Sm
            Mass Number = 147
            Relative Atomic Mass = 146.9149044(19)
            Isotopic Composition = 0.1499(18)
            Standard Atomic Weight = 150.36(2)
            Notes = g

            Atomic Number = 62
            Atomic Symbol = Sm
            Mass Number = 148
            Relative Atomic Mass = 147.9148292(19)
            Isotopic Composition = 0.1124(10)
            Standard Atomic Weight = 150.36(2)
            Notes = g

            Atomic Number = 62
            Atomic Symbol = Sm
            Mass Number = 149
            Relative Atomic Mass = 148.9171921(18)
            Isotopic Composition = 0.1382(7)
            Standard Atomic Weight = 150.36(2)
            Notes = g

            Atomic Number = 62
            Atomic Symbol = Sm
            Mass Number = 150
            Relative Atomic Mass = 149.9172829(18)
            Isotopic Composition = 0.0738(1)
            Standard Atomic Weight = 150.36(2)
            Notes = g

            Atomic Number = 62
            Atomic Symbol = Sm
            Mass Number = 152
            Relative Atomic Mass = 151.9197397(18)
            Isotopic Composition = 0.2675(16)
            Standard Atomic Weight = 150.36(2)
            Notes = g

            Atomic Number = 62
            Atomic Symbol = Sm
            Mass Number = 154
            Relative Atomic Mass = 153.9222169(20)
            Isotopic Composition = 0.2275(29)
            Standard Atomic Weight = 150.36(2)
            Notes = g

            Atomic Number = 63
            Atomic Symbol = Eu
            Mass Number = 151
            Relative Atomic Mass = 150.9198578(18)
            Isotopic Composition = 0.4781(6)
            Standard Atomic Weight = 151.964(1)
            Notes = g

            Atomic Number = 63
            Atomic Symbol = Eu
            Mass Number = 153
            Relative Atomic Mass = 152.9212380(18)
            Isotopic Composition = 0.5219(6)
            Standard Atomic Weight = 151.964(1)
            Notes = g

            Atomic Number = 64
            Atomic Symbol = Gd
            Mass Number = 152
            Relative Atomic Mass = 151.9197995(18)
            Isotopic Composition = 0.0020(1)
            Standard Atomic Weight = 157.25(3)
            Notes = g

            Atomic Number = 64
            Atomic Symbol = Gd
            Mass Number = 154
            Relative Atomic Mass = 153.9208741(17)
            Isotopic Composition = 0.0218(3)
            Standard Atomic Weight = 157.25(3)
            Notes = g

            Atomic Number = 64
            Atomic Symbol = Gd
            Mass Number = 155
            Relative Atomic Mass = 154.9226305(17)
            Isotopic Composition = 0.1480(12)
            Standard Atomic Weight = 157.25(3)
            Notes = g

            Atomic Number = 64
            Atomic Symbol = Gd
            Mass Number = 156
            Relative Atomic Mass = 155.9221312(17)
            Isotopic Composition = 0.2047(9)
            Standard Atomic Weight = 157.25(3)
            Notes = g

            Atomic Number = 64
            Atomic Symbol = Gd
            Mass Number = 157
            Relative Atomic Mass = 156.9239686(17)
            Isotopic Composition = 0.1565(2)
            Standard Atomic Weight = 157.25(3)
            Notes = g

            Atomic Number = 64
            Atomic Symbol = Gd
            Mass Number = 158
            Relative Atomic Mass = 157.9241123(17)
            Isotopic Composition = 0.2484(7)
            Standard Atomic Weight = 157.25(3)
            Notes = g

            Atomic Number = 64
            Atomic Symbol = Gd
            Mass Number = 160
            Relative Atomic Mass = 159.9270624(18)
            Isotopic Composition = 0.2186(19)
            Standard Atomic Weight = 157.25(3)
            Notes = g

            Atomic Number = 65
            Atomic Symbol = Tb
            Mass Number = 159
            Relative Atomic Mass = 158.9253547(19)
            Isotopic Composition = 1
            Standard Atomic Weight = 158.92535(2)
            Notes = &nbsp;

            Atomic Number = 66
            Atomic Symbol = Dy
            Mass Number = 156
            Relative Atomic Mass = 155.9242847(17)
            Isotopic Composition = 0.00056(3)
            Standard Atomic Weight = 162.500(1)
            Notes = g

            Atomic Number = 66
            Atomic Symbol = Dy
            Mass Number = 158
            Relative Atomic Mass = 157.9244159(31)
            Isotopic Composition = 0.00095(3)
            Standard Atomic Weight = 162.500(1)
            Notes = g

            Atomic Number = 66
            Atomic Symbol = Dy
            Mass Number = 160
            Relative Atomic Mass = 159.9252046(20)
            Isotopic Composition = 0.02329(18)
            Standard Atomic Weight = 162.500(1)
            Notes = g

            Atomic Number = 66
            Atomic Symbol = Dy
            Mass Number = 161
            Relative Atomic Mass = 160.9269405(20)
            Isotopic Composition = 0.18889(42)
            Standard Atomic Weight = 162.500(1)
            Notes = g

            Atomic Number = 66
            Atomic Symbol = Dy
            Mass Number = 162
            Relative Atomic Mass = 161.9268056(20)
            Isotopic Composition = 0.25475(36)
            Standard Atomic Weight = 162.500(1)
            Notes = g

            Atomic Number = 66
            Atomic Symbol = Dy
            Mass Number = 163
            Relative Atomic Mass = 162.9287383(20)
            Isotopic Composition = 0.24896(42)
            Standard Atomic Weight = 162.500(1)
            Notes = g

            Atomic Number = 66
            Atomic Symbol = Dy
            Mass Number = 164
            Relative Atomic Mass = 163.9291819(20)
            Isotopic Composition = 0.28260(54)
            Standard Atomic Weight = 162.500(1)
            Notes = g

            Atomic Number = 67
            Atomic Symbol = Ho
            Mass Number = 165
            Relative Atomic Mass = 164.9303288(21)
            Isotopic Composition = 1
            Standard Atomic Weight = 164.93033(2)
            Notes = &nbsp;

            Atomic Number = 68
            Atomic Symbol = Er
            Mass Number = 162
            Relative Atomic Mass = 161.9287884(20)
            Isotopic Composition = 0.00139(5)
            Standard Atomic Weight = 167.259(3)
            Notes = g

            Atomic Number = 68
            Atomic Symbol = Er
            Mass Number = 164
            Relative Atomic Mass = 163.9292088(20)
            Isotopic Composition = 0.01601(3)
            Standard Atomic Weight = 167.259(3)
            Notes = g

            Atomic Number = 68
            Atomic Symbol = Er
            Mass Number = 166
            Relative Atomic Mass = 165.9302995(22)
            Isotopic Composition = 0.33503(36)
            Standard Atomic Weight = 167.259(3)
            Notes = g

            Atomic Number = 68
            Atomic Symbol = Er
            Mass Number = 167
            Relative Atomic Mass = 166.9320546(22)
            Isotopic Composition = 0.22869(9)
            Standard Atomic Weight = 167.259(3)
            Notes = g

            Atomic Number = 68
            Atomic Symbol = Er
            Mass Number = 168
            Relative Atomic Mass = 167.9323767(22)
            Isotopic Composition = 0.26978(18)
            Standard Atomic Weight = 167.259(3)
            Notes = g

            Atomic Number = 68
            Atomic Symbol = Er
            Mass Number = 170
            Relative Atomic Mass = 169.9354702(26)
            Isotopic Composition = 0.14910(36)
            Standard Atomic Weight = 167.259(3)
            Notes = g

            Atomic Number = 69
            Atomic Symbol = Tm
            Mass Number = 169
            Relative Atomic Mass = 168.9342179(22)
            Isotopic Composition = 1
            Standard Atomic Weight = 168.93422(2)
            Notes = &nbsp;

            Atomic Number = 70
            Atomic Symbol = Yb
            Mass Number = 168
            Relative Atomic Mass = 167.9338896(22)
            Isotopic Composition = 0.00123(3)
            Standard Atomic Weight = 173.054(5)
            Notes = g

            Atomic Number = 70
            Atomic Symbol = Yb
            Mass Number = 170
            Relative Atomic Mass = 169.9347664(22)
            Isotopic Composition = 0.02982(39)
            Standard Atomic Weight = 173.054(5)
            Notes = g

            Atomic Number = 70
            Atomic Symbol = Yb
            Mass Number = 171
            Relative Atomic Mass = 170.9363302(22)
            Isotopic Composition = 0.1409(14)
            Standard Atomic Weight = 173.054(5)
            Notes = g

            Atomic Number = 70
            Atomic Symbol = Yb
            Mass Number = 172
            Relative Atomic Mass = 171.9363859(22)
            Isotopic Composition = 0.2168(13)
            Standard Atomic Weight = 173.054(5)
            Notes = g

            Atomic Number = 70
            Atomic Symbol = Yb
            Mass Number = 173
            Relative Atomic Mass = 172.9382151(22)
            Isotopic Composition = 0.16103(63)
            Standard Atomic Weight = 173.054(5)
            Notes = g

            Atomic Number = 70
            Atomic Symbol = Yb
            Mass Number = 174
            Relative Atomic Mass = 173.9388664(22)
            Isotopic Composition = 0.32026(80)
            Standard Atomic Weight = 173.054(5)
            Notes = g

            Atomic Number = 70
            Atomic Symbol = Yb
            Mass Number = 176
            Relative Atomic Mass = 175.9425764(24)
            Isotopic Composition = 0.12996(83)
            Standard Atomic Weight = 173.054(5)
            Notes = g

            Atomic Number = 71
            Atomic Symbol = Lu
            Mass Number = 175
            Relative Atomic Mass = 174.9407752(20)
            Isotopic Composition = 0.97401(13)
            Standard Atomic Weight = 174.9668(1)
            Notes = g

            Atomic Number = 71
            Atomic Symbol = Lu
            Mass Number = 176
            Relative Atomic Mass = 175.9426897(20)
            Isotopic Composition = 0.02599(13)
            Standard Atomic Weight = 174.9668(1)
            Notes = g

            Atomic Number = 72
            Atomic Symbol = Hf
            Mass Number = 174
            Relative Atomic Mass = 173.9400461(28)
            Isotopic Composition = 0.0016(1)
            Standard Atomic Weight = 178.49(2)
            Notes = &nbsp;

            Atomic Number = 72
            Atomic Symbol = Hf
            Mass Number = 176
            Relative Atomic Mass = 175.9414076(22)
            Isotopic Composition = 0.0526(7)
            Standard Atomic Weight = 178.49(2)
            Notes = &nbsp;

            Atomic Number = 72
            Atomic Symbol = Hf
            Mass Number = 177
            Relative Atomic Mass = 176.9432277(20)
            Isotopic Composition = 0.1860(9)
            Standard Atomic Weight = 178.49(2)
            Notes = &nbsp;

            Atomic Number = 72
            Atomic Symbol = Hf
            Mass Number = 178
            Relative Atomic Mass = 177.9437058(20)
            Isotopic Composition = 0.2728(7)
            Standard Atomic Weight = 178.49(2)
            Notes = &nbsp;

            Atomic Number = 72
            Atomic Symbol = Hf
            Mass Number = 179
            Relative Atomic Mass = 178.9458232(20)
            Isotopic Composition = 0.1362(2)
            Standard Atomic Weight = 178.49(2)
            Notes = &nbsp;

            Atomic Number = 72
            Atomic Symbol = Hf
            Mass Number = 180
            Relative Atomic Mass = 179.9465570(20)
            Isotopic Composition = 0.3508(16)
            Standard Atomic Weight = 178.49(2)
            Notes = &nbsp;

            Atomic Number = 73
            Atomic Symbol = Ta
            Mass Number = 180
            Relative Atomic Mass = 179.9474648(24)
            Isotopic Composition = 0.0001201(32)
            Standard Atomic Weight = 180.94788(2)
            Notes = &nbsp;

            Atomic Number = 73
            Atomic Symbol = Ta
            Mass Number = 181
            Relative Atomic Mass = 180.9479958(20)
            Isotopic Composition = 0.9998799(32)
            Standard Atomic Weight = 180.94788(2)
            Notes = &nbsp;

            Atomic Number = 74
            Atomic Symbol = W
            Mass Number = 180
            Relative Atomic Mass = 179.9467108(20)
            Isotopic Composition = 0.0012(1)
            Standard Atomic Weight = 183.84(1)
            Notes = &nbsp;

            Atomic Number = 74
            Atomic Symbol = W
            Mass Number = 182
            Relative Atomic Mass = 181.94820394(91)
            Isotopic Composition = 0.2650(16)
            Standard Atomic Weight = 183.84(1)
            Notes = &nbsp;

            Atomic Number = 74
            Atomic Symbol = W
            Mass Number = 183
            Relative Atomic Mass = 182.95022275(90)
            Isotopic Composition = 0.1431(4)
            Standard Atomic Weight = 183.84(1)
            Notes = &nbsp;

            Atomic Number = 74
            Atomic Symbol = W
            Mass Number = 184
            Relative Atomic Mass = 183.95093092(94)
            Isotopic Composition = 0.3064(2)
            Standard Atomic Weight = 183.84(1)
            Notes = &nbsp;

            Atomic Number = 74
            Atomic Symbol = W
            Mass Number = 186
            Relative Atomic Mass = 185.9543628(17)
            Isotopic Composition = 0.2843(19)
            Standard Atomic Weight = 183.84(1)
            Notes = &nbsp;

            Atomic Number = 75
            Atomic Symbol = Re
            Mass Number = 185
            Relative Atomic Mass = 184.9529545(13)
            Isotopic Composition = 0.3740(2)
            Standard Atomic Weight = 186.207(1)
            Notes = &nbsp;

            Atomic Number = 75
            Atomic Symbol = Re
            Mass Number = 187
            Relative Atomic Mass = 186.9557501(16)
            Isotopic Composition = 0.6260(2)
            Standard Atomic Weight = 186.207(1)
            Notes = &nbsp;

            Atomic Number = 76
            Atomic Symbol = Os
            Mass Number = 184
            Relative Atomic Mass = 183.9524885(14)
            Isotopic Composition = 0.0002(1)
            Standard Atomic Weight = 190.23(3)
            Notes = g

            Atomic Number = 76
            Atomic Symbol = Os
            Mass Number = 186
            Relative Atomic Mass = 185.9538350(16)
            Isotopic Composition = 0.0159(3)
            Standard Atomic Weight = 190.23(3)
            Notes = g

            Atomic Number = 76
            Atomic Symbol = Os
            Mass Number = 187
            Relative Atomic Mass = 186.9557474(16)
            Isotopic Composition = 0.0196(2)
            Standard Atomic Weight = 190.23(3)
            Notes = g

            Atomic Number = 76
            Atomic Symbol = Os
            Mass Number = 188
            Relative Atomic Mass = 187.9558352(16)
            Isotopic Composition = 0.1324(8)
            Standard Atomic Weight = 190.23(3)
            Notes = g

            Atomic Number = 76
            Atomic Symbol = Os
            Mass Number = 189
            Relative Atomic Mass = 188.9581442(17)
            Isotopic Composition = 0.1615(5)
            Standard Atomic Weight = 190.23(3)
            Notes = g

            Atomic Number = 76
            Atomic Symbol = Os
            Mass Number = 190
            Relative Atomic Mass = 189.9584437(17)
            Isotopic Composition = 0.2626(2)
            Standard Atomic Weight = 190.23(3)
            Notes = g

            Atomic Number = 76
            Atomic Symbol = Os
            Mass Number = 192
            Relative Atomic Mass = 191.9614770(29)
            Isotopic Composition = 0.4078(19)
            Standard Atomic Weight = 190.23(3)
            Notes = g

            Atomic Number = 77
            Atomic Symbol = Ir
            Mass Number = 191
            Relative Atomic Mass = 190.9605893(21)
            Isotopic Composition = 0.373(2)
            Standard Atomic Weight = 192.217(3)
            Notes = &nbsp;

            Atomic Number = 77
            Atomic Symbol = Ir
            Mass Number = 193
            Relative Atomic Mass = 192.9629216(21)
            Isotopic Composition = 0.627(2)
            Standard Atomic Weight = 192.217(3)
            Notes = &nbsp;

            Atomic Number = 78
            Atomic Symbol = Pt
            Mass Number = 190
            Relative Atomic Mass = 189.9599297(63)
            Isotopic Composition = 0.00012(2)
            Standard Atomic Weight = 195.084(9)
            Notes = &nbsp;

            Atomic Number = 78
            Atomic Symbol = Pt
            Mass Number = 192
            Relative Atomic Mass = 191.9610387(32)
            Isotopic Composition = 0.00782(24)
            Standard Atomic Weight = 195.084(9)
            Notes = &nbsp;

            Atomic Number = 78
            Atomic Symbol = Pt
            Mass Number = 194
            Relative Atomic Mass = 193.9626809(10)
            Isotopic Composition = 0.3286(40)
            Standard Atomic Weight = 195.084(9)
            Notes = &nbsp;

            Atomic Number = 78
            Atomic Symbol = Pt
            Mass Number = 195
            Relative Atomic Mass = 194.9647917(10)
            Isotopic Composition = 0.3378(24)
            Standard Atomic Weight = 195.084(9)
            Notes = &nbsp;

            Atomic Number = 78
            Atomic Symbol = Pt
            Mass Number = 196
            Relative Atomic Mass = 195.96495209(99)
            Isotopic Composition = 0.2521(34)
            Standard Atomic Weight = 195.084(9)
            Notes = &nbsp;

            Atomic Number = 78
            Atomic Symbol = Pt
            Mass Number = 198
            Relative Atomic Mass = 197.9678949(23)
            Isotopic Composition = 0.07356(130)
            Standard Atomic Weight = 195.084(9)
            Notes = &nbsp;

            Atomic Number = 79
            Atomic Symbol = Au
            Mass Number = 197
            Relative Atomic Mass = 196.96656879(71)
            Isotopic Composition = 1
            Standard Atomic Weight = 196.966569(5)
            Notes = &nbsp;

            Atomic Number = 80
            Atomic Symbol = Hg
            Mass Number = 196
            Relative Atomic Mass = 195.9658326(32)
            Isotopic Composition = 0.0015(1)
            Standard Atomic Weight = 200.592(3)
            Notes = &nbsp;

            Atomic Number = 80
            Atomic Symbol = Hg
            Mass Number = 198
            Relative Atomic Mass = 197.96676860(52)
            Isotopic Composition = 0.0997(20)
            Standard Atomic Weight = 200.592(3)
            Notes = &nbsp;

            Atomic Number = 80
            Atomic Symbol = Hg
            Mass Number = 199
            Relative Atomic Mass = 198.96828064(46)
            Isotopic Composition = 0.1687(22)
            Standard Atomic Weight = 200.592(3)
            Notes = &nbsp;

            Atomic Number = 80
            Atomic Symbol = Hg
            Mass Number = 200
            Relative Atomic Mass = 199.96832659(47)
            Isotopic Composition = 0.2310(19)
            Standard Atomic Weight = 200.592(3)
            Notes = &nbsp;

            Atomic Number = 80
            Atomic Symbol = Hg
            Mass Number = 201
            Relative Atomic Mass = 200.97030284(69)
            Isotopic Composition = 0.1318(9)
            Standard Atomic Weight = 200.592(3)
            Notes = &nbsp;

            Atomic Number = 80
            Atomic Symbol = Hg
            Mass Number = 202
            Relative Atomic Mass = 201.97064340(69)
            Isotopic Composition = 0.2986(26)
            Standard Atomic Weight = 200.592(3)
            Notes = &nbsp;

            Atomic Number = 80
            Atomic Symbol = Hg
            Mass Number = 204
            Relative Atomic Mass = 203.97349398(53)
            Isotopic Composition = 0.0687(15)
            Standard Atomic Weight = 200.592(3)
            Notes = &nbsp;

            Atomic Number = 81
            Atomic Symbol = Tl
            Mass Number = 203
            Relative Atomic Mass = 202.9723446(14)
            Isotopic Composition = 0.2952(1)
            Standard Atomic Weight = [204.382,204.385]
            Notes = &nbsp;

            Atomic Number = 81
            Atomic Symbol = Tl
            Mass Number = 205
            Relative Atomic Mass = 204.9744278(14)
            Isotopic Composition = 0.7048(1)
            Standard Atomic Weight = [204.382,204.385]
            Notes = &nbsp;

            Atomic Number = 82
            Atomic Symbol = Pb
            Mass Number = 204
            Relative Atomic Mass = 203.9730440(13)
            Isotopic Composition = 0.014(1)
            Standard Atomic Weight = 207.2(1)
            Notes = g,r

            Atomic Number = 82
            Atomic Symbol = Pb
            Mass Number = 206
            Relative Atomic Mass = 205.9744657(13)
            Isotopic Composition = 0.241(1)
            Standard Atomic Weight = 207.2(1)
            Notes = g,r

            Atomic Number = 82
            Atomic Symbol = Pb
            Mass Number = 207
            Relative Atomic Mass = 206.9758973(13)
            Isotopic Composition = 0.221(1)
            Standard Atomic Weight = 207.2(1)
            Notes = g,r

            Atomic Number = 82
            Atomic Symbol = Pb
            Mass Number = 208
            Relative Atomic Mass = 207.9766525(13)
            Isotopic Composition = 0.524(1)
            Standard Atomic Weight = 207.2(1)
            Notes = g,r

            Atomic Number = 83
            Atomic Symbol = Bi
            Mass Number = 209
            Relative Atomic Mass = 208.9803991(16)
            Isotopic Composition = 1
            Standard Atomic Weight = 208.98040(1)
            Notes = &nbsp;

            Atomic Number = 84
            Atomic Symbol = Po
            Mass Number = 209
            Relative Atomic Mass = 208.9824308(20)
            Isotopic Composition =
            Standard Atomic Weight = [209]
            Notes = &nbsp;

            Atomic Number = 84
            Atomic Symbol = Po
            Mass Number = 210
            Relative Atomic Mass = 209.9828741(13)
            Isotopic Composition =
            Standard Atomic Weight = [209]
            Notes = &nbsp;

            Atomic Number = 85
            Atomic Symbol = At
            Mass Number = 210
            Relative Atomic Mass = 209.9871479(83)
            Isotopic Composition =
            Standard Atomic Weight = [210]
            Notes = &nbsp;

            Atomic Number = 85
            Atomic Symbol = At
            Mass Number = 211
            Relative Atomic Mass = 210.9874966(30)
            Isotopic Composition =
            Standard Atomic Weight = [210]
            Notes = &nbsp;

            Atomic Number = 86
            Atomic Symbol = Rn
            Mass Number = 211
            Relative Atomic Mass = 210.9906011(73)
            Isotopic Composition =
            Standard Atomic Weight = [222]
            Notes = &nbsp;

            Atomic Number = 86
            Atomic Symbol = Rn
            Mass Number = 220
            Relative Atomic Mass = 220.0113941(23)
            Isotopic Composition =
            Standard Atomic Weight = [222]
            Notes = &nbsp;

            Atomic Number = 86
            Atomic Symbol = Rn
            Mass Number = 222
            Relative Atomic Mass = 222.0175782(25)
            Isotopic Composition =
            Standard Atomic Weight = [222]
            Notes = &nbsp;

            Atomic Number = 87
            Atomic Symbol = Fr
            Mass Number = 223
            Relative Atomic Mass = 223.0197360(25)
            Isotopic Composition =
            Standard Atomic Weight = [223]
            Notes = &nbsp;

            Atomic Number = 88
            Atomic Symbol = Ra
            Mass Number = 223
            Relative Atomic Mass = 223.0185023(27)
            Isotopic Composition =
            Standard Atomic Weight = [226]
            Notes = &nbsp;

            Atomic Number = 88
            Atomic Symbol = Ra
            Mass Number = 224
            Relative Atomic Mass = 224.0202120(23)
            Isotopic Composition =
            Standard Atomic Weight = [226]
            Notes = &nbsp;

            Atomic Number = 88
            Atomic Symbol = Ra
            Mass Number = 226
            Relative Atomic Mass = 226.0254103(25)
            Isotopic Composition =
            Standard Atomic Weight = [226]
            Notes = &nbsp;

            Atomic Number = 88
            Atomic Symbol = Ra
            Mass Number = 228
            Relative Atomic Mass = 228.0310707(26)
            Isotopic Composition =
            Standard Atomic Weight = [226]
            Notes = &nbsp;

            Atomic Number = 89
            Atomic Symbol = Ac
            Mass Number = 227
            Relative Atomic Mass = 227.0277523(25)
            Isotopic Composition =
            Standard Atomic Weight = [227]
            Notes = &nbsp;

            Atomic Number = 90
            Atomic Symbol = Th
            Mass Number = 230
            Relative Atomic Mass = 230.0331341(19)
            Isotopic Composition =
            Standard Atomic Weight = 232.0377(4)
            Notes = g

            Atomic Number = 90
            Atomic Symbol = Th
            Mass Number = 232
            Relative Atomic Mass = 232.0380558(21)
            Isotopic Composition = 1
            Standard Atomic Weight = 232.0377(4)
            Notes = g

            Atomic Number = 91
            Atomic Symbol = Pa
            Mass Number = 231
            Relative Atomic Mass = 231.0358842(24)
            Isotopic Composition = 1
            Standard Atomic Weight = 231.03588(2)
            Notes = &nbsp;

            Atomic Number = 92
            Atomic Symbol = U
            Mass Number = 233
            Relative Atomic Mass = 233.0396355(29)
            Isotopic Composition =
            Standard Atomic Weight = 238.02891(3)
            Notes = g,m

            Atomic Number = 92
            Atomic Symbol = U
            Mass Number = 234
            Relative Atomic Mass = 234.0409523(19)
            Isotopic Composition = 0.000054(5)
            Standard Atomic Weight = 238.02891(3)
            Notes = g,m

            Atomic Number = 92
            Atomic Symbol = U
            Mass Number = 235
            Relative Atomic Mass = 235.0439301(19)
            Isotopic Composition = 0.007204(6)
            Standard Atomic Weight = 238.02891(3)
            Notes = g,m

            Atomic Number = 92
            Atomic Symbol = U
            Mass Number = 236
            Relative Atomic Mass = 236.0455682(19)
            Isotopic Composition =
            Standard Atomic Weight = 238.02891(3)
            Notes = g,m

            Atomic Number = 92
            Atomic Symbol = U
            Mass Number = 238
            Relative Atomic Mass = 238.0507884(20)
            Isotopic Composition = 0.992742(10)
            Standard Atomic Weight = 238.02891(3)
            Notes = g,m

            Atomic Number = 93
            Atomic Symbol = Np
            Mass Number = 236
            Relative Atomic Mass = 236.046570(54)
            Isotopic Composition =
            Standard Atomic Weight = [237]
            Notes = &nbsp;

            Atomic Number = 93
            Atomic Symbol = Np
            Mass Number = 237
            Relative Atomic Mass = 237.0481736(19)
            Isotopic Composition =
            Standard Atomic Weight = [237]
            Notes = &nbsp;

            Atomic Number = 94
            Atomic Symbol = Pu
            Mass Number = 238
            Relative Atomic Mass = 238.0495601(19)
            Isotopic Composition =
            Standard Atomic Weight = [244]
            Notes = &nbsp;

            Atomic Number = 94
            Atomic Symbol = Pu
            Mass Number = 239
            Relative Atomic Mass = 239.0521636(19)
            Isotopic Composition =
            Standard Atomic Weight = [244]
            Notes = &nbsp;

            Atomic Number = 94
            Atomic Symbol = Pu
            Mass Number = 240
            Relative Atomic Mass = 240.0538138(19)
            Isotopic Composition =
            Standard Atomic Weight = [244]
            Notes = &nbsp;

            Atomic Number = 94
            Atomic Symbol = Pu
            Mass Number = 241
            Relative Atomic Mass = 241.0568517(19)
            Isotopic Composition =
            Standard Atomic Weight = [244]
            Notes = &nbsp;

            Atomic Number = 94
            Atomic Symbol = Pu
            Mass Number = 242
            Relative Atomic Mass = 242.0587428(20)
            Isotopic Composition =
            Standard Atomic Weight = [244]
            Notes = &nbsp;

            Atomic Number = 94
            Atomic Symbol = Pu
            Mass Number = 244
            Relative Atomic Mass = 244.0642053(56)
            Isotopic Composition =
            Standard Atomic Weight = [244]
            Notes = &nbsp;

            Atomic Number = 95
            Atomic Symbol = Am
            Mass Number = 241
            Relative Atomic Mass = 241.0568293(19)
            Isotopic Composition =
            Standard Atomic Weight = &nbsp;
            Notes = &nbsp;

            Atomic Number = 95
            Atomic Symbol = Am
            Mass Number = 243
            Relative Atomic Mass = 243.0613813(24)
            Isotopic Composition =
            Standard Atomic Weight = &nbsp;
            Notes = &nbsp;

            Atomic Number = 96
            Atomic Symbol = Cm
            Mass Number = 243
            Relative Atomic Mass = 243.0613893(22)
            Isotopic Composition =
            Standard Atomic Weight = &nbsp;
            Notes = &nbsp;

            Atomic Number = 96
            Atomic Symbol = Cm
            Mass Number = 244
            Relative Atomic Mass = 244.0627528(19)
            Isotopic Composition =
            Standard Atomic Weight = &nbsp;
            Notes = &nbsp;

            Atomic Number = 96
            Atomic Symbol = Cm
            Mass Number = 245
            Relative Atomic Mass = 245.0654915(22)
            Isotopic Composition =
            Standard Atomic Weight = &nbsp;
            Notes = &nbsp;

            Atomic Number = 96
            Atomic Symbol = Cm
            Mass Number = 246
            Relative Atomic Mass = 246.0672238(22)
            Isotopic Composition =
            Standard Atomic Weight = &nbsp;
            Notes = &nbsp;

            Atomic Number = 96
            Atomic Symbol = Cm
            Mass Number = 247
            Relative Atomic Mass = 247.0703541(47)
            Isotopic Composition =
            Standard Atomic Weight = &nbsp;
            Notes = &nbsp;

            Atomic Number = 96
            Atomic Symbol = Cm
            Mass Number = 248
            Relative Atomic Mass = 248.0723499(56)
            Isotopic Composition =
            Standard Atomic Weight = &nbsp;
            Notes = &nbsp;

            Atomic Number = 97
            Atomic Symbol = Bk
            Mass Number = 247
            Relative Atomic Mass = 247.0703073(59)
            Isotopic Composition =
            Standard Atomic Weight = &nbsp;
            Notes = &nbsp;

            Atomic Number = 97
            Atomic Symbol = Bk
            Mass Number = 249
            Relative Atomic Mass = 249.0749877(27)
            Isotopic Composition =
            Standard Atomic Weight = &nbsp;
            Notes = &nbsp;

            Atomic Number = 98
            Atomic Symbol = Cf
            Mass Number = 249
            Relative Atomic Mass = 249.0748539(23)
            Isotopic Composition =
            Standard Atomic Weight = &nbsp;
            Notes = &nbsp;

            Atomic Number = 98
            Atomic Symbol = Cf
            Mass Number = 250
            Relative Atomic Mass = 250.0764062(22)
            Isotopic Composition =
            Standard Atomic Weight = &nbsp;
            Notes = &nbsp;

            Atomic Number = 98
            Atomic Symbol = Cf
            Mass Number = 251
            Relative Atomic Mass = 251.0795886(48)
            Isotopic Composition =
            Standard Atomic Weight = &nbsp;
            Notes = &nbsp;

            Atomic Number = 98
            Atomic Symbol = Cf
            Mass Number = 252
            Relative Atomic Mass = 252.0816272(56)
            Isotopic Composition =
            Standard Atomic Weight = &nbsp;
            Notes = &nbsp;

            Atomic Number = 99
            Atomic Symbol = Es
            Mass Number = 252
            Relative Atomic Mass = 252.082980(54)
            Isotopic Composition =
            Standard Atomic Weight = &nbsp;
            Notes = &nbsp;

            Atomic Number = 100
            Atomic Symbol = Fm
            Mass Number = 257
            Relative Atomic Mass = 257.0951061(69)
            Isotopic Composition =
            Standard Atomic Weight = &nbsp;
            Notes = &nbsp;

            Atomic Number = 101
            Atomic Symbol = Md
            Mass Number = 258
            Relative Atomic Mass = 258.0984315(50)
            Isotopic Composition =
            Standard Atomic Weight = &nbsp;
            Notes = &nbsp;

            Atomic Number = 101
            Atomic Symbol = Md
            Mass Number = 260
            Relative Atomic Mass = 260.10365(34#)
            Isotopic Composition =
            Standard Atomic Weight = &nbsp;
            Notes = &nbsp;

            Atomic Number = 102
            Atomic Symbol = No
            Mass Number = 259
            Relative Atomic Mass = 259.10103(11#)
            Isotopic Composition =
            Standard Atomic Weight = &nbsp;
            Notes = &nbsp;

            Atomic Number = 103
            Atomic Symbol = Lr
            Mass Number = 262
            Relative Atomic Mass = 262.10961(22#)
            Isotopic Composition =
            Standard Atomic Weight = &nbsp;
            Notes = &nbsp;

            Atomic Number = 104
            Atomic Symbol = Rf
            Mass Number = 267
            Relative Atomic Mass = 267.12179(62#)
            Isotopic Composition =
            Standard Atomic Weight = &nbsp;
            Notes = &nbsp;

            Atomic Number = 105
            Atomic Symbol = Db
            Mass Number = 268
            Relative Atomic Mass = 268.12567(57#)
            Isotopic Composition =
            Standard Atomic Weight = &nbsp;
            Notes = &nbsp;

            Atomic Number = 106
            Atomic Symbol = Sg
            Mass Number = 271
            Relative Atomic Mass = 271.13393(63#)
            Isotopic Composition =
            Standard Atomic Weight = &nbsp;
            Notes = &nbsp;

            Atomic Number = 107
            Atomic Symbol = Bh
            Mass Number = 272
            Relative Atomic Mass = 272.13826(58#)
            Isotopic Composition =
            Standard Atomic Weight = &nbsp;
            Notes = &nbsp;

            Atomic Number = 108
            Atomic Symbol = Hs
            Mass Number = 270
            Relative Atomic Mass = 270.13429(27#)
            Isotopic Composition =
            Standard Atomic Weight = &nbsp;
            Notes = &nbsp;

            Atomic Number = 109
            Atomic Symbol = Mt
            Mass Number = 276
            Relative Atomic Mass = 276.15159(59#)
            Isotopic Composition =
            Standard Atomic Weight = &nbsp;
            Notes = &nbsp;

            Atomic Number = 110
            Atomic Symbol = Ds
            Mass Number = 281
            Relative Atomic Mass = 281.16451(59#)
            Isotopic Composition =
            Standard Atomic Weight = &nbsp;
            Notes = &nbsp;

            Atomic Number = 111
            Atomic Symbol = Rg
            Mass Number = 280
            Relative Atomic Mass = 280.16514(61#)
            Isotopic Composition =
            Standard Atomic Weight = &nbsp;
            Notes = &nbsp;

            Atomic Number = 112
            Atomic Symbol = Cn
            Mass Number = 285
            Relative Atomic Mass = 285.17712(60#)
            Isotopic Composition =
            Standard Atomic Weight = &nbsp;
            Notes = &nbsp;

            Atomic Number = 113
            Atomic Symbol = Nh
            Mass Number = 284
            Relative Atomic Mass = 284.17873(62#)
            Isotopic Composition =
            Standard Atomic Weight = &nbsp;
            Notes = &nbsp;

            Atomic Number = 114
            Atomic Symbol = Fl
            Mass Number = 289
            Relative Atomic Mass = 289.19042(60#)
            Isotopic Composition =
            Standard Atomic Weight = &nbsp;
            Notes = &nbsp;

            Atomic Number = 115
            Atomic Symbol = Mc
            Mass Number = 288
            Relative Atomic Mass = 288.19274(62#)
            Isotopic Composition =
            Standard Atomic Weight = &nbsp;
            Notes = &nbsp;

            Atomic Number = 116
            Atomic Symbol = Lv
            Mass Number = 293
            Relative Atomic Mass = 293.20449(60#)
            Isotopic Composition =
            Standard Atomic Weight = &nbsp;
            Notes = &nbsp;

            Atomic Number = 117
            Atomic Symbol = Ts
            Mass Number = 292
            Relative Atomic Mass = 292.20746(75#)
            Isotopic Composition =
            Standard Atomic Weight = &nbsp;
            Notes = &nbsp;

            Atomic Number = 118
            Atomic Symbol = Og
            Mass Number = 294
            Relative Atomic Mass = 294.21392(71#)
            Isotopic Composition =
            Standard Atomic Weight = &nbsp;
            Notes = &nbsp;";
    }
}