/**
 * Theme Initialization Script
 * ===========================
 * 
 * This script handles early theme initialization to prevent FOUC (Flash of Unstyled Content)
 * and FOUT (Flash of Incorrect Theme). It must be loaded as early as possible in the document head.
 * 
 * Features:
 * - Applies saved theme preference from localStorage
 * - Falls back to system preference via prefers-color-scheme
 * - Handles errors gracefully with proper fallbacks
 * - Sets appropriate classes and data attributes for CSS targeting
 * 
 * Author: SOAPY Documentation System
 * Updated: 2025-01-25
 */

(function() {
  try {
    // Apply the theme as quickly as possible to avoid flash
    let saved;
    try {
      saved = localStorage.getItem('theme');
    } catch (storageError) {
      // Handle localStorage access errors (private browsing, cookies disabled)
      console.warn('Unable to access localStorage:', storageError);
    }

    let prefersDark = false;
    try {
      prefersDark = window.matchMedia && window.matchMedia('(prefers-color-scheme: dark)').matches;
    } catch (mediaError) {
      // Handle matchMedia access errors (unsupported browsers)
      console.warn('Unable to check media query:', mediaError);
    }

    const theme = saved ? saved : (prefersDark ? 'dark' : 'light');
    document.documentElement.setAttribute('data-theme', theme);
    // Add theme class for CSS targeting
    document.documentElement.classList.add(theme + '-theme');
    // Mark when theme is set for CSS transitions
    document.documentElement.classList.add('theme-set');
  } catch (error) {
    // Fallback for any other unexpected errors
    console.warn('Theme initialization failed:', error);
    // Apply default theme as fallback
    try {
      document.documentElement.setAttribute('data-theme', 'light');
      document.documentElement.classList.add('light-theme');
      document.documentElement.classList.add('theme-set');
    } catch (fallbackError) {
      // Silent catch for extreme cases
      console.error('Critical theme fallback failed:', fallbackError);
    }
  }
})();