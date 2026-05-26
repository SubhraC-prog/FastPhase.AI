"""
Gradient Program Optimizer Module

This module provides advanced gradient optimization algorithms for HPLC method development.
It implements linear solvent strength theory, design space exploration, and multi-parameter
optimization for robust chromatographic separations.

Author: Chromatography AI Team
Version: 1.0.0
References: Snyder, Dolan, Neue, Horváth, Schoenmakers
"""

import logging
import numpy as np
from typing import Dict, List, Tuple, Optional, Any, Callable
from dataclasses import dataclass, field
from enum import Enum
import math
from io import BytesIO
import base64
from datetime import datetime
import json

try:
    from scipy.optimize import minimize, differential_evolution
    from scipy.interpolate import interp1d
    _HAS_SCIPY = True
except ImportError:
    _HAS_SCIPY = False

# matplotlib is imported lazily inside plot methods to avoid crashes
# in headless/server environments (no display available).
_plt = None

def _get_plt():
    """Lazy matplotlib import to avoid headless crash."""
    global _plt
    if _plt is None:
        try:
            import matplotlib
            matplotlib.use('Agg')  # non-interactive backend safe for servers
            import matplotlib.pyplot as plt
            _plt = plt
        except ImportError:
            pass
    return _plt

# Configure logging
logger = logging.getLogger(__name__)


class GradientType(Enum):
    """Enum for gradient types"""
    LINEAR = "linear"
    STEP = "step"
    MULTI_LINEAR = "multi_linear"
    EXPONENTIAL = "exponential"
    CONVEX = "convex"
    CONCAVE = "concave"


class OptimizationObjective(Enum):
    """Enum for optimization objectives"""
    RESOLUTION = "resolution"
    ANALYSIS_TIME = "analysis_time"
    PEAK_CAPACITY = "peak_capacity"
    ROBUSTNESS = "robustness"
    BALANCED = "balanced"


@dataclass
class GradientSegment:
    """Data class for a single gradient segment"""
    start_time: float
    end_time: float
    start_b: float
    end_b: float
    type: GradientType = GradientType.LINEAR
    curve_factor: float = 1.0  # For non-linear gradients


@dataclass
class GradientProgram:
    """Data class for complete gradient program"""
    segments: List[GradientSegment]
    flow_rate: float
    column_temperature: float
    equilibration_time: float
    injection_volume: float
    total_runtime: float
    
    def to_table(self) -> str:
        """Convert to formatted gradient table"""
        table = "Time (min)\t%B\tEvent\n"
        for segment in self.segments:
            table += f"{segment.start_time:.2f}\t{segment.start_b:.1f}\tStart segment\n"
            table += f"{segment.end_time:.2f}\t{segment.end_b:.1f}\tEnd segment\n"
        return table
    
    def to_dict(self) -> Dict:
        """Convert to dictionary"""
        return {
            'segments': [{
                'start_time': s.start_time,
                'end_time': s.end_time,
                'start_b': s.start_b,
                'end_b': s.end_b,
                'type': s.type.value,
                'curve_factor': s.curve_factor
            } for s in self.segments],
            'flow_rate': self.flow_rate,
            'column_temperature': self.column_temperature,
            'equilibration_time': self.equilibration_time,
            'injection_volume': self.injection_volume,
            'total_runtime': self.total_runtime
        }


@dataclass
class PeakPrediction:
    """Data class for predicted peak parameters"""
    compound_name: str
    retention_time: float
    peak_width: float
    asymmetry: float
    efficiency: float
    resolution_from_previous: Optional[float] = None


@dataclass
class OptimizationResult:
    """Data class for optimization results"""
    gradient_program: GradientProgram
    peak_predictions: List[PeakPrediction]
    objective_value: float
    constraints_satisfied: bool
    design_space: Dict[str, Tuple[float, float]]
    robustness_metrics: Dict[str, float]
    confidence_score: float
    warnings: List[str] = field(default_factory=list)
    
    def summary(self) -> Dict[str, Any]:
        """Generate summary dictionary"""
        return {
            'total_runtime': self.gradient_program.total_runtime,
            'peak_capacity': self._calculate_peak_capacity(),
            'min_resolution': min([p.resolution_from_previous for p in self.peak_predictions[1:]] 
                                  if len(self.peak_predictions) > 1 else 1.5),
            'objective_value': self.objective_value,
            'confidence': self.confidence_score,
            'warnings': self.warnings
        }
    
    def _calculate_peak_capacity(self) -> float:
        """Calculate peak capacity for the gradient"""
        if len(self.peak_predictions) < 2:
            return 0
        first_rt = self.peak_predictions[0].retention_time
        last_rt = self.peak_predictions[-1].retention_time
        avg_width = np.mean([p.peak_width for p in self.peak_predictions])
        return (last_rt - first_rt) / avg_width


class GradientOptimizer:
    """
    Advanced gradient program optimizer for HPLC method development.
    
    This class implements various optimization algorithms based on:
    - Linear Solvent Strength Theory (Snyder, Dolan)
    - Design Space Exploration (Quality by Design)
    - Multi-parameter optimization
    - Robustness analysis
    
    References:
    [1] Snyder LR, Dolan JW. High-Performance Gradient Elution. Wiley; 2007.
    [2] Neue UD. HPLC Columns: Theory, Technology, and Practice. Wiley-VCH; 1997.
    [3] Schoenmakers PJ. Optimization of Chromatographic Selectivity. Elsevier; 1986.
    [4] Horváth C, Melander W. Theory of chromatography. J Chromatogr Sci. 1977;15(9):393-404.
    [5] Dolan JW. Temperature selectivity in reversed-phase high performance liquid chromatography.
        J Chromatogr A. 2002;965(1-2):195-205.
    """
    
    # Literature-derived constants
    S_VALUE_RANGE = (3, 8)  # Slope of log k' vs φ
    K0_RANGE = (1, 20)  # Retention factor at 0% B
    MIN_PRESSURE = 50  # bar
    MAX_PRESSURE = 400  # bar (standard HPLC)
    OPTIMAL_K_RANGE = (1, 20)  # Optimal retention factor range
    
    def __init__(self, 
                 column_length: float = 150,  # mm
                 column_id: float = 4.6,  # mm
                 particle_size: float = 3.5,  # µm
                 dead_time: Optional[float] = None,
                 system_volume: float = 1.2,  # mL (gradient delay volume)
                 max_pressure: float = 400):  # bar
        """
        Initialize the gradient optimizer with column parameters.
        
        Args:
            column_length: Column length in mm
            column_id: Column internal diameter in mm
            particle_size: Particle size in µm
            dead_time: Column dead time (t0) in minutes (calculated if None)
            system_volume: Gradient delay volume in mL
            max_pressure: Maximum system pressure in bar
        """
        self.column_length = column_length
        self.column_id = column_id
        self.particle_size = particle_size
        self.system_volume = system_volume
        self.max_pressure = max_pressure
        
        # Calculate column volume and dead time
        self.column_volume = self._calculate_column_volume()
        
        # Set dead time based on typical flow rate
        if dead_time is None:
            self.dead_time = self.column_volume / 2.0  # Assuming 2 mL/min flow
        else:
            self.dead_time = dead_time
        
        logger.info(f"Initialized optimizer with column volume: {self.column_volume:.2f} mL")
    
    def _calculate_column_volume(self) -> float:
        """Calculate geometric column volume in mL"""
        radius_cm = (self.column_id / 2) / 10  # Convert mm to cm
        length_cm = self.column_length / 10  # Convert mm to cm
        volume_mL = math.pi * radius_cm**2 * length_cm
        return volume_mL
    
    def optimize_gradient(self,
                         compounds: List[Dict[str, Any]],
                         objective: OptimizationObjective = OptimizationObjective.BALANCED,
                         gradient_type: GradientType = GradientType.LINEAR,
                         constraints: Optional[Dict] = None) -> OptimizationResult:
        """
        Optimize gradient program for given compounds.
        
        Args:
            compounds: List of compound dictionaries with properties:
                - name: Compound identifier
                - logp: LogP value
                - pka: pKa value (if ionizable)
                - charge: Formal charge at pH
                - hydrophobicity: Hydrophobic parameter
            objective: Optimization objective
            gradient_type: Type of gradient to optimize
            constraints: Additional constraints (max time, min resolution, etc.)
            
        Returns:
            OptimizationResult containing optimized gradient and predictions
        """
        logger.info(f"Starting gradient optimization for {len(compounds)} compounds")
        
        # Set default constraints
        if constraints is None:
            constraints = {
                'max_time': 60,
                'min_resolution': 1.5,
                'max_pressure': self.max_pressure,
                'min_b': 5,
                'max_b': 95
            }
        
        # Estimate S-values (slope of log k' vs φ) for each compound
        # Based on molecular weight and hydrophobicity [Snyder, 2007]
        for compound in compounds:
            compound['S'] = self._estimate_s_value(compound)
            compound['logk0'] = self._estimate_logk0(compound)
        
        # Perform optimization based on gradient type
        if gradient_type == GradientType.LINEAR:
            result = self._optimize_linear_gradient(compounds, objective, constraints)
        elif gradient_type == GradientType.MULTI_LINEAR:
            result = self._optimize_multilinear_gradient(compounds, objective, constraints)
        else:
            result = self._optimize_general_gradient(compounds, objective, constraints, gradient_type)
        
        # Calculate robustness metrics
        result.robustness_metrics = self._calculate_robustness(result.gradient_program, compounds)
        
        # Calculate confidence score
        result.confidence_score = self._calculate_confidence(result)
        
        logger.info(f"Optimization complete. Runtime: {result.gradient_program.total_runtime:.1f} min")
        
        return result
    
    def _estimate_s_value(self, compound: Dict) -> float:
        """
        Estimate S-value (slope of log k' vs φ) based on compound properties.
        
        Based on Snyder-Dolan relationship: S ≈ 0.48 * MW^0.44 [Snyder, 2007]
        """
        mw = compound.get('molecular_weight', 300)
        logp = compound.get('logp', 3)
        
        # Base S from molecular weight
        S = 0.48 * (mw ** 0.44)
        
        # Adjust for hydrophobicity
        S = S * (1 + 0.1 * (logp - 3))
        
        # Adjust for ionizable groups
        if compound.get('charge', 0) != 0:
            S = S * 0.9  # Ionizable compounds show smaller S
        
        return np.clip(S, self.S_VALUE_RANGE[0], self.S_VALUE_RANGE[1])
    
    def _estimate_logk0(self, compound: Dict) -> float:
        """
        Estimate log k'0 (retention factor at 0% B).
        
        Based on relationship with logP [Horváth, 1977]
        """
        logp = compound.get('logp', 3)
        
        # Base logk0 from logP
        logk0 = 0.5 * logp
        
        # Adjust for polar groups
        tpsa = compound.get('tpsa', 0)
        logk0 = logk0 - 0.01 * tpsa
        
        return np.clip(logk0, -1, 3)
    
    def _optimize_linear_gradient(self, 
                                  compounds: List[Dict],
                                  objective: OptimizationObjective,
                                  constraints: Dict) -> OptimizationResult:
        """
        Optimize simple linear gradient (initial B, final B, time).
        
        Based on Linear Solvent Strength Theory [Snyder, 2007]
        """
        
        def objective_function(params):
            """Objective function for optimization"""
            initial_b, final_b, gradient_time = params
            
            # Clip parameters to valid ranges
            initial_b = np.clip(initial_b, constraints['min_b'], constraints['max_b'] - 10)
            final_b = np.clip(final_b, initial_b + 10, constraints['max_b'])
            gradient_time = np.clip(gradient_time, 5, constraints['max_time'])
            
            # Create gradient program
            gradient = self._create_linear_gradient(initial_b, final_b, gradient_time)
            
            # Predict peaks
            peaks = self._predict_peaks(gradient, compounds)
            
            # Calculate objective value
            obj_val = self._calculate_objective(peaks, objective, constraints)
            
            return obj_val
        
        # Initial guess based on compound properties
        avg_logp = np.mean([c.get('logp', 3) for c in compounds])
        initial_guess = [
            self._estimate_initial_b(avg_logp),
            self._estimate_final_b(avg_logp),
            self._estimate_gradient_time(compounds)
        ]
        
        # Bounds for parameters
        bounds = [
            (constraints['min_b'], constraints['max_b'] - 10),
            (constraints['min_b'] + 10, constraints['max_b']),
            (5, constraints['max_time'])
        ]
        
        # Perform optimization
        result = minimize(
            objective_function,
            initial_guess,
            bounds=bounds,
            method='L-BFGS-B'
        )
        
        # Create final gradient
        initial_b, final_b, gradient_time = result.x
        gradient = self._create_linear_gradient(initial_b, final_b, gradient_time)
        
        # Predict peaks
        peaks = self._predict_peaks(gradient, compounds)
        
        return OptimizationResult(
            gradient_program=gradient,
            peak_predictions=peaks,
            objective_value=result.fun,
            constraints_satisfied=self._check_constraints(peaks, constraints),
            design_space=self._calculate_design_space(compounds),
            robustness_metrics={},
            confidence_score=0.0
        )
    
    def _optimize_multilinear_gradient(self,
                                      compounds: List[Dict],
                                      objective: OptimizationObjective,
                                      constraints: Dict) -> OptimizationResult:
        """
        Optimize multi-linear gradient with multiple segments.
        
        Implements design space exploration [Neue, 1997]
        """
        
        def objective_function(params):
            """Objective function with multiple segments"""
            n_segments = 3  # Optimize for 3 segments
            segment_params = params.reshape(n_segments, 3)  # Each: [start_b, end_b, time]
            
            # Build gradient segments
            segments = []
            current_time = 0
            for i, (start_b, end_b, time) in enumerate(segment_params):
                start_b = np.clip(start_b, constraints['min_b'], constraints['max_b'])
                end_b = np.clip(end_b, start_b, constraints['max_b'])
                time = np.clip(time, 2, constraints['max_time'] / n_segments)
                
                segments.append(GradientSegment(
                    start_time=current_time,
                    end_time=current_time + time,
                    start_b=start_b,
                    end_b=end_b
                ))
                current_time += time
            
            # Create gradient program
            gradient = GradientProgram(
                segments=segments,
                flow_rate=2.0,  # Default flow
                column_temperature=30,  # Default temp
                equilibration_time=5,
                injection_volume=5,
                total_runtime=current_time + 5
            )
            
            # Predict peaks
            peaks = self._predict_peaks(gradient, compounds)
            
            return self._calculate_objective(peaks, objective, constraints)
        
        # Initial guess based on compound distribution
        n_params = 9  # 3 segments × 3 parameters
        initial_guess = np.array([30, 50, 10, 50, 70, 10, 70, 95, 10])  # Example
        
        bounds = [(constraints['min_b'], constraints['max_b'])] * 6 + [(2, 20)] * 3
        
        # Use differential evolution for global optimization
        result = differential_evolution(
            objective_function,
            bounds,
            maxiter=100,
            popsize=15,
            seed=42
        )
        
        # Reconstruct gradient from optimized parameters
        segments = []
        current_time = 0
        opt_params = result.x.reshape(3, 3)
        
        for i, (start_b, end_b, time) in enumerate(opt_params):
            segments.append(GradientSegment(
                start_time=current_time,
                end_time=current_time + time,
                start_b=start_b,
                end_b=end_b
            ))
            current_time += time
        
        gradient = GradientProgram(
            segments=segments,
            flow_rate=2.0,
            column_temperature=30,
            equilibration_time=5,
            injection_volume=5,
            total_runtime=current_time + 5
        )
        
        peaks = self._predict_peaks(gradient, compounds)
        
        return OptimizationResult(
            gradient_program=gradient,
            peak_predictions=peaks,
            objective_value=result.fun,
            constraints_satisfied=self._check_constraints(peaks, constraints),
            design_space=self._calculate_design_space(compounds),
            robustness_metrics={},
            confidence_score=0.0
        )
    
    def _optimize_general_gradient(self,
                                  compounds: List[Dict],
                                  objective: OptimizationObjective,
                                  constraints: Dict,
                                  gradient_type: GradientType) -> OptimizationResult:
        """
        Optimize general non-linear gradient shapes.
        
        Based on theory of non-linear gradients [Schoenmakers, 1986]
        """
        # Implementation for exponential, convex, concave gradients
        # This is a simplified version - full implementation would include
        # proper modeling of non-linear solvent strength
        
        def objective_function(params):
            if gradient_type == GradientType.EXPONENTIAL:
                initial_b, final_b, time, exponent = params
            else:
                initial_b, final_b, time, curve_factor = params
            
            # Create non-linear gradient
            gradient = self._create_nonlinear_gradient(
                initial_b, final_b, time, gradient_type, 
                curve_factor if gradient_type != GradientType.EXPONENTIAL else exponent
            )
            
            peaks = self._predict_peaks(gradient, compounds)
            return self._calculate_objective(peaks, objective, constraints)
        
        # Simplified optimization for non-linear gradients
        avg_logp = np.mean([c.get('logp', 3) for c in compounds])
        initial_guess = [
            self._estimate_initial_b(avg_logp),
            self._estimate_final_b(avg_logp),
            self._estimate_gradient_time(compounds),
            1.5  # Curve factor
        ]
        
        bounds = [
            (constraints['min_b'], constraints['max_b'] - 10),
            (constraints['min_b'] + 10, constraints['max_b']),
            (5, constraints['max_time']),
            (0.5, 3.0)  # Curve factor range
        ]
        
        result = minimize(objective_function, initial_guess, bounds=bounds, method='L-BFGS-B')
        
        # Create gradient
        gradient = self._create_nonlinear_gradient(
            result.x[0], result.x[1], result.x[2], gradient_type, result.x[3]
        )
        
        peaks = self._predict_peaks(gradient, compounds)
        
        return OptimizationResult(
            gradient_program=gradient,
            peak_predictions=peaks,
            objective_value=result.fun,
            constraints_satisfied=self._check_constraints(peaks, constraints),
            design_space=self._calculate_design_space(compounds),
            robustness_metrics={},
            confidence_score=0.0
        )
    
    def _create_linear_gradient(self, initial_b: float, final_b: float, 
                                gradient_time: float) -> GradientProgram:
        """Create a simple linear gradient program"""
        segment = GradientSegment(
            start_time=0,
            end_time=gradient_time,
            start_b=initial_b,
            end_b=final_b,
            type=GradientType.LINEAR
        )
        
        return GradientProgram(
            segments=[segment],
            flow_rate=2.0,  # mL/min (typical)
            column_temperature=30,  # °C
            equilibration_time=5,  # min
            injection_volume=5,  # µL
            total_runtime=gradient_time + 5
        )
    
    def _create_nonlinear_gradient(self, initial_b: float, final_b: float,
                                   gradient_time: float, gradient_type: GradientType,
                                   curve_factor: float) -> GradientProgram:
        """Create a non-linear gradient program"""
        segment = GradientSegment(
            start_time=0,
            end_time=gradient_time,
            start_b=initial_b,
            end_b=final_b,
            type=gradient_type,
            curve_factor=curve_factor
        )
        
        return GradientProgram(
            segments=[segment],
            flow_rate=2.0,
            column_temperature=30,
            equilibration_time=5,
            injection_volume=5,
            total_runtime=gradient_time + 5
        )
    
    def _predict_peaks(self, gradient: GradientProgram, 
                       compounds: List[Dict]) -> List[PeakPrediction]:
        """
        Predict retention times and peak shapes using LSS theory.
        
        Implements Snyder's linear solvent strength theory [Snyder, 2007]
        and Neue's gradient elution equations [Neue, 1997]
        """
        predictions = []
        
        for i, compound in enumerate(compounds):
            # Calculate retention time using gradient equation
            # tR = t0 + (t0/b) * log(1 + b * k0) where b = S * Δφ * t0 / tG
            t0 = self.dead_time
            
            # Get compound parameters
            S = compound.get('S', 5)
            logk0 = compound.get('logk0', 1)
            k0 = 10 ** logk0
            
            # Calculate gradient steepness parameter (b)
            delta_phi = (gradient.segments[-1].end_b - gradient.segments[0].start_b) / 100
            tG = gradient.total_runtime - gradient.equilibration_time
            
            b = S * delta_phi * t0 / tG
            
            # Calculate retention time
            if b > 0:
                tR = t0 + (t0 / b) * math.log(1 + b * k0)
            else:
                tR = t0 * (1 + k0)
            
            # Account for system volume (gradient delay)
            tR = tR + self.system_volume / gradient.flow_rate
            
            # Calculate peak width [Neue, 1997]
            N = self._estimate_plate_count(compound)
            peak_width = tR / math.sqrt(N) * 4  # 4σ width
            
            # Calculate asymmetry (based on compound properties)
            asymmetry = self._estimate_asymmetry(compound)
            
            # Calculate efficiency
            efficiency = N * (1 + k0) / (1 + k0 * t0 / tR)  # Gradient efficiency
            
            prediction = PeakPrediction(
                compound_name=compound.get('name', f'Compound_{i+1}'),
                retention_time=tR,
                peak_width=peak_width,
                asymmetry=asymmetry,
                efficiency=efficiency
            )
            
            predictions.append(prediction)
        
        # Sort by retention time
        predictions.sort(key=lambda x: x.retention_time)
        
        # Calculate resolutions between adjacent peaks
        for i in range(1, len(predictions)):
            prev = predictions[i-1]
            curr = predictions[i]
            
            resolution = (curr.retention_time - prev.retention_time) / \
                        ((curr.peak_width + prev.peak_width) / 2)
            curr.resolution_from_previous = resolution
        
        return predictions
    
    def _estimate_plate_count(self, compound: Dict) -> float:
        """
        Estimate theoretical plate count for compound.
        
        Based on column parameters and compound properties [Snyder, 2007]
        """
        # Base plate count from column
        N_base = 300 * (self.column_length / self.particle_size)
        
        # Adjust for compound retention
        logp = compound.get('logp', 3)
        N = N_base * (1 + 0.1 * (logp - 3))
        
        return N
    
    def _estimate_asymmetry(self, compound: Dict) -> float:
        """
        Estimate peak asymmetry factor.
        
        Based on compound's interaction with stationary phase
        """
        base_asymmetry = 1.0
        
        # Adjust for basic compounds (silanol interactions)
        if compound.get('charge', 0) > 0:
            base_asymmetry += 0.3
        
        # Adjust for very hydrophobic compounds
        if compound.get('logp', 3) > 5:
            base_asymmetry += 0.2
        
        return base_asymmetry
    
    def _estimate_initial_b(self, avg_logp: float) -> float:
        """Estimate initial %B based on average logP"""
        if avg_logp < 1:
            return 10
        elif avg_logp < 3:
            return 25
        elif avg_logp < 5:
            return 40
        else:
            return 55
    
    def _estimate_final_b(self, avg_logp: float) -> float:
        """Estimate final %B based on average logP"""
        initial = self._estimate_initial_b(avg_logp)
        return min(95, initial + 40)
    
    def _estimate_gradient_time(self, compounds: List[Dict]) -> float:
        """Estimate gradient time based on compound diversity"""
        # Calculate range of hydrophobicities
        logp_values = [c.get('logp', 3) for c in compounds]
        logp_range = max(logp_values) - min(logp_values)
        
        # Base time on hydrophobicity range
        base_time = 20 + logp_range * 5
        
        # Adjust for number of compounds
        base_time = base_time + len(compounds) * 2
        
        return min(60, base_time)
    
    def _calculate_objective(self, peaks: List[PeakPrediction],
                            objective: OptimizationObjective,
                            constraints: Dict) -> float:
        """
        Calculate objective function value for optimization.
        
        Implements multi-objective optimization [Schoenmakers, 1986]
        """
        if objective == OptimizationObjective.RESOLUTION:
            if len(peaks) < 2:
                return -10  # Penalty for single peak
            resolutions = [p.resolution_from_previous for p in peaks[1:]]
            min_res = min(resolutions) if resolutions else 1.5
            return -min_res  # Negative for minimization
        
        elif objective == OptimizationObjective.ANALYSIS_TIME:
            last_peak = peaks[-1] if peaks else None
            if last_peak:
                return last_peak.retention_time
            return 60
        
        elif objective == OptimizationObjective.PEAK_CAPACITY:
            if len(peaks) < 2:
                return -10
            first_rt = peaks[0].retention_time
            last_rt = peaks[-1].retention_time
            avg_width = np.mean([p.peak_width for p in peaks])
            peak_capacity = (last_rt - first_rt) / avg_width
            return -peak_capacity
        
        elif objective == OptimizationObjective.ROBUSTNESS:
            # Maximize robustness (minimize sensitivity to parameter changes)
            return -self._calculate_robustness_score(peaks)
        
        else:  # BALANCED
            # Weighted combination of objectives
            if len(peaks) < 2:
                return 100
            
            # Resolution component
            resolutions = [p.resolution_from_previous for p in peaks[1:]]
            min_res = min(resolutions) if resolutions else 1.5
            res_score = max(0, 2.0 - min_res) * 10  # Penalize below 2.0
            
            # Time component
            last_rt = peaks[-1].retention_time
            time_score = last_rt / 10  # Normalize
            
            return res_score + time_score
    
    def _calculate_robustness_score(self, peaks: List[PeakPrediction]) -> float:
        """Calculate robustness score based on peak spacing"""
        if len(peaks) < 2:
            return 1.0
        
        resolutions = [p.resolution_from_previous for p in peaks[1:]]
        
        # Calculate standard deviation of resolutions
        res_std = np.std(resolutions)
        
        # Lower std indicates more robust method
        return 1.0 / (1.0 + res_std)
    
    def _check_constraints(self, peaks: List[PeakPrediction],
                          constraints: Dict) -> bool:
        """Check if constraints are satisfied"""
        if not peaks:
            return False
        
        # Check minimum resolution
        if len(peaks) > 1:
            resolutions = [p.resolution_from_previous for p in peaks[1:]]
            if min(resolutions) < constraints.get('min_resolution', 1.5):
                return False
        
        # Check maximum time
        last_rt = peaks[-1].retention_time
        if last_rt > constraints.get('max_time', 60):
            return False
        
        return True
    
    def _calculate_design_space(self, compounds: List[Dict]) -> Dict[str, Tuple[float, float]]:
        """
        Calculate design space for robust operation.
        
        Implements Quality by Design principles [USP <1220>]
        """
        design_space = {
            'ph_range': (2.0, 8.0),  # Typical range for silica-based columns
            'temperature_range': (20, 40),  # °C
            'flow_rate_range': (1.0, 3.0),  # mL/min
            'composition_tolerance': (95, 105)  # % of optimal
        }
        
        # Adjust based on compound properties
        acidic_pka = [c.get('pka_acidic', 14) for c in compounds if 'pka_acidic' in c]
        basic_pka = [c.get('pka_basic', 0) for c in compounds if 'pka_basic' in c]
        
        if acidic_pka:
            # Stay at least 1 unit below lowest pKa
            design_space['ph_range'] = (2.0, min(acidic_pka) - 1)
        
        if basic_pka:
            # Stay at least 1 unit above highest pKa
            design_space['ph_range'] = (max(basic_pka) + 1, 8.0)
        
        return design_space
    
    def _calculate_robustness(self, gradient: GradientProgram,
                            compounds: List[Dict]) -> Dict[str, float]:
        """
        Calculate robustness metrics by Monte Carlo simulation.
        
        Evaluates method robustness to small parameter variations [USP <621>]
        """
        n_simulations = 100
        base_peaks = self._predict_peaks(gradient, compounds)
        base_resolutions = [p.resolution_from_previous for p in base_peaks[1:]] if len(base_peaks) > 1 else [2.0]
        
        # Parameters to vary
        variations = {
            'flow_rate': 0.1,  # 10% variation
            'temperature': 2,  # ±2°C
            'initial_b': 1,  # ±1%
            'final_b': 1,  # ±1%
        }
        
        resolution_variations = []
        retention_time_variations = []
        
        for _ in range(n_simulations):
            # Create perturbed gradient
            perturbed_gradient = self._perturb_gradient(gradient, variations)
            
            # Predict peaks
            peaks = self._predict_peaks(perturbed_gradient, compounds)
            
            if len(peaks) > 1:
                resolutions = [p.resolution_from_previous for p in peaks[1:]]
                resolution_variations.append(np.std(resolutions))
            
            retention_times = [p.retention_time for p in peaks]
            retention_time_variations.append(np.std(retention_times))
        
        return {
            'resolution_robustness': 1.0 - np.mean(resolution_variations) if resolution_variations else 0,
            'retention_time_robustness': 1.0 / (1.0 + np.mean(retention_time_variations)) if retention_time_variations else 0,
            'critical_pair_robustness': self._calculate_critical_pair_robustness(gradient, compounds)
        }
    
    def _perturb_gradient(self, gradient: GradientProgram,
                         variations: Dict) -> GradientProgram:
        """Apply random perturbations to gradient parameters"""
        import random
        
        # Copy segments
        new_segments = []
        for segment in gradient.segments:
            new_segment = GradientSegment(
                start_time=segment.start_time,
                end_time=segment.end_time,
                start_b=segment.start_b + random.uniform(-variations['initial_b'], variations['initial_b']),
                end_b=segment.end_b + random.uniform(-variations['final_b'], variations['final_b']),
                type=segment.type,
                curve_factor=segment.curve_factor
            )
            new_segments.append(new_segment)
        
        return GradientProgram(
            segments=new_segments,
            flow_rate=gradient.flow_rate * random.uniform(0.9, 1.1),
            column_temperature=gradient.column_temperature + random.uniform(-2, 2),
            equilibration_time=gradient.equilibration_time,
            injection_volume=gradient.injection_volume,
            total_runtime=gradient.total_runtime
        )
    
    def _calculate_critical_pair_robustness(self, gradient: GradientProgram,
                                          compounds: List[Dict]) -> float:
        """Calculate robustness specifically for critical peak pairs"""
        peaks = self._predict_peaks(gradient, compounds)
        
        if len(peaks) < 2:
            return 1.0
        
        # Find critical pair (closest peaks)
        min_resolution = float('inf')
        critical_idx = 0
        
        for i in range(1, len(peaks)):
            if peaks[i].resolution_from_previous and peaks[i].resolution_from_previous < min_resolution:
                min_resolution = peaks[i].resolution_from_previous
                critical_idx = i
        
        # Analyze sensitivity of critical pair
        sensitivities = []
        
        # Test flow rate variation
        for flow_mult in [0.95, 1.05]:
            test_gradient = GradientProgram(
                segments=gradient.segments,
                flow_rate=gradient.flow_rate * flow_mult,
                column_temperature=gradient.column_temperature,
                equilibration_time=gradient.equilibration_time,
                injection_volume=gradient.injection_volume,
                total_runtime=gradient.total_runtime
            )
            test_peaks = self._predict_peaks(test_gradient, compounds)
            if len(test_peaks) > critical_idx:
                sensitivities.append(abs(test_peaks[critical_idx].resolution_from_previous - min_resolution))
        
        return 1.0 / (1.0 + np.mean(sensitivities)) if sensitivities else 1.0
    
    def _calculate_confidence(self, result: OptimizationResult) -> float:
        """Calculate overall confidence score for the optimization"""
        confidence = 0.9  # Base confidence
        
        # Reduce confidence based on potential issues
        if len(result.peak_predictions) < 2:
            confidence *= 0.8
        
        # Check resolution
        resolutions = [p.resolution_from_previous for p in result.peak_predictions[1:]] if len(result.peak_predictions) > 1 else [2.0]
        if min(resolutions) < 1.2:
            confidence *= 0.7
        elif min(resolutions) < 1.5:
            confidence *= 0.9
        
        # Check retention times
        last_rt = result.peak_predictions[-1].retention_time if result.peak_predictions else 0
        if last_rt > result.gradient_program.total_runtime - 2:
            confidence *= 0.9  # Risk of eluting after gradient
        
        # Consider robustness metrics
        if result.robustness_metrics:
            avg_robustness = np.mean(list(result.robustness_metrics.values()))
            confidence *= avg_robustness
        
        return round(confidence, 3)
    
    def generate_gradient_table(self, gradient: GradientProgram) -> str:
        """Generate formatted gradient table for method submission"""
        table = []
        table.append("=" * 60)
        table.append("GRADIENT PROGRAM")
        table.append("=" * 60)
        table.append(f"Column: {self.column_length} × {self.column_id} mm, {self.particle_size} µm")
        table.append(f"Flow Rate: {gradient.flow_rate:.1f} mL/min")
        table.append(f"Temperature: {gradient.column_temperature}°C")
        table.append(f"Injection Volume: {gradient.injection_volume} µL")
        table.append(f"Equilibration Time: {gradient.equilibration_time} min")
        table.append("-" * 60)
        table.append("Time (min)\t%B\tEvent")
        table.append("-" * 60)
        
        for segment in gradient.segments:
            table.append(f"{segment.start_time:.2f}\t\t{segment.start_b:.1f}\tStart")
            table.append(f"{segment.end_time:.2f}\t\t{segment.end_b:.1f}\tEnd")
        
        table.append("-" * 60)
        table.append(f"Total Runtime: {gradient.total_runtime:.1f} min")
        table.append("=" * 60)
        
        return "\n".join(table)
    
    def plot_gradient_profile(self, gradient: GradientProgram) -> str:
        """Generate base64 encoded plot of gradient profile"""
        try:
            plt = _get_plt()
            if plt is None:
                logger.warning("matplotlib not available; skipping gradient plot.")
                return ""

            plt.figure(figsize=(10, 6))
            
            # Generate time points
            t_points = []
            b_points = []
            
            for segment in gradient.segments:
                t_seg = np.linspace(segment.start_time, segment.end_time, 50)
                
                if segment.type == GradientType.LINEAR:
                    b_seg = np.interp(t_seg, [segment.start_time, segment.end_time], 
                                     [segment.start_b, segment.end_b])
                elif segment.type == GradientType.EXPONENTIAL:
                    b_seg = segment.start_b + (segment.end_b - segment.start_b) * \
                           (1 - np.exp(-segment.curve_factor * (t_seg - segment.start_time) / 
                                      (segment.end_time - segment.start_time)))
                else:
                    b_seg = np.interp(t_seg, [segment.start_time, segment.end_time], 
                                     [segment.start_b, segment.end_b])
                
                t_points.extend(t_seg)
                b_points.extend(b_seg)
            
            plt.plot(t_points, b_points, 'b-', linewidth=2)
            plt.xlabel('Time (min)')
            plt.ylabel('%B')
            plt.title('Gradient Profile')
            plt.grid(True, alpha=0.3)
            plt.ylim(0, 100)
            
            # Add markers at segment boundaries
            for segment in gradient.segments:
                plt.plot(segment.start_time, segment.start_b, 'ro', markersize=8)
                plt.plot(segment.end_time, segment.end_b, 'go', markersize=8)
            
            # Convert plot to base64 string
            buffer = BytesIO()
            plt.savefig(buffer, format='png', dpi=100, bbox_inches='tight')
            buffer.seek(0)
            image_base64 = base64.b64encode(buffer.getvalue()).decode()
            plt.close()
            
            return image_base64
            
        except Exception as e:
            logger.warning(f"Could not generate plot: {e}")
            return ""


# Utility functions for common gradient calculations

def calculate_gradient_steepness(initial_b: float, final_b: float, 
                                gradient_time: float) -> float:
    """Calculate gradient steepness (%B/min)"""
    return (final_b - initial_b) / gradient_time


def calculate_peak_capacity(gradient_time: float, dead_time: float,
                          plate_count: float) -> float:
    """
    Calculate theoretical peak capacity [Neue, 1997]
    
    Pc = 1 + (tG / (t0 * 4)) * sqrt(N)
    """
    return 1 + (gradient_time / (dead_time * 4)) * math.sqrt(plate_count)


def calculate_resolution(k1: float, k2: float, N: float) -> float:
    """
    Calculate resolution between two peaks [Snyder, 2007]
    
    Rs = (sqrt(N)/4) * ((α-1)/α) * (k2/(1+k2))
    """
    alpha = k2 / k1
    k_avg = (k1 + k2) / 2
    
    return (math.sqrt(N) / 4) * ((alpha - 1) / alpha) * (k_avg / (1 + k_avg))


def estimate_required_gradient_time(logp_range: float, 
                                   n_compounds: int,
                                   desired_resolution: float = 1.5) -> float:
    """Estimate required gradient time based on compound diversity"""
    base_time = 20 + logp_range * 5
    compound_factor = 1 + (n_compounds - 2) * 0.2  # 20% increase per extra compound
    resolution_factor = desired_resolution / 1.5
    
    return base_time * compound_factor * resolution_factor


# Example usage and test
if __name__ == "__main__":
    # Test compounds (pharmaceutical mixture)
    test_compounds = [
        {
            'name': 'Paracetamol',
            'logp': 0.46,
            'molecular_weight': 151,
            'tpsa': 49.3,
            'pka_acidic': 9.5,
            'pka_basic': 1.7
        },
        {
            'name': 'Caffeine',
            'logp': -0.07,
            'molecular_weight': 194,
            'tpsa': 58.4,
            'pka_basic': 10.4
        },
        {
            'name': 'Aspirin',
            'logp': 1.19,
            'molecular_weight': 180,
            'tpsa': 63.6,
            'pka_acidic': 3.5
        },
        {
            'name': 'Ibuprofen',
            'logp': 3.97,
            'molecular_weight': 206,
            'tpsa': 37.3,
            'pka_acidic': 4.4
        }
    ]
    
    print("=" * 70)
    print("GRADIENT OPTIMIZER TEST")
    print("=" * 70)
    
    # Initialize optimizer
    optimizer = GradientOptimizer(
        column_length=150,
        column_id=4.6,
        particle_size=3.5,
        system_volume=1.2
    )
    
    print(f"\nColumn volume: {optimizer.column_volume:.2f} mL")
    print(f"Dead time: {optimizer.dead_time:.2f} min")
    
    # Test different optimization objectives
    objectives = [
        OptimizationObjective.BALANCED,
        OptimizationObjective.RESOLUTION,
        OptimizationObjective.ANALYSIS_TIME
    ]
    
    for obj in objectives:
        print(f"\n{'-' * 50}")
        print(f"Optimizing for: {obj.value}")
        print(f"{'-' * 50}")
        
        result = optimizer.optimize_gradient(
            compounds=test_compounds,
            objective=obj,
            gradient_type=GradientType.LINEAR
        )
        
        summary = result.summary()
        print(f"Runtime: {summary['total_runtime']:.1f} min")
        print(f"Peak capacity: {summary['peak_capacity']:.1f}")
        print(f"Min resolution: {summary['min_resolution']:.2f}")
        print(f"Confidence: {summary['confidence']:.2f}")
        
        # Print gradient table
        print("\n" + optimizer.generate_gradient_table(result.gradient_program))
    
    # Test multi-linear gradient
    print(f"\n{'-' * 50}")
    print("Multi-linear gradient optimization")
    print(f"{'-' * 50}")
    
    result = optimizer.optimize_gradient(
        compounds=test_compounds,
        objective=OptimizationObjective.BALANCED,
        gradient_type=GradientType.MULTI_LINEAR
    )
    
    print(f"Runtime: {result.gradient_program.total_runtime:.1f} min")
    print("Segments:")
    for i, seg in enumerate(result.gradient_program.segments, 1):
        print(f"  Segment {i}: {seg.start_b:.1f}% → {seg.end_b:.1f}% over {seg.end_time - seg.start_time:.1f} min")
    
    print("\n" + "=" * 70)
    print("Test complete")
    print("=" * 70)