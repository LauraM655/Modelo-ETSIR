
<script>
document.addEventListener("DOMContentLoaded", function() {
  
  // 🔹 1. INTERACTIVIDAD DEL ÍNDICE ACORDEÓN (PRINCIPAL Y SUBTEMAS)
  function setupAccordion(tabContainer) {
    const indexLinks = tabContainer.querySelectorAll('.index-box a');
    const contentCards = tabContainer.querySelectorAll('.content-card');
    
    if (contentCards.length > 0) {
      contentCards[0].classList.add('active');
      if (indexLinks[0]) indexLinks[0].classList.add('active');
    }
    
    indexLinks.forEach(link => {
      link.addEventListener('click', function(e) {
        e.preventDefault();
        const targetId = this.getAttribute('href').substring(1);
        const targetSection = document.getElementById(targetId);
        
        if (!targetSection) return;
        
        contentCards.forEach(card => card.classList.remove('active'));
        indexLinks.forEach(l => l.classList.remove('active'));
        
        targetSection.classList.add('active');
        this.classList.add('active');
        
        targetSection.scrollIntoView({ behavior: 'smooth', block: 'start' });
      });
    });
  }
  
  const tabPanes = document.querySelectorAll('.tab-pane, .section.level2');
  tabPanes.forEach(pane => {
    if (pane.querySelector('.index-box')) {
      setupAccordion(pane);
    }
  });
  
  const mainTabs = document.querySelectorAll('.nav-tabs a[data-toggle="tab"]');
  mainTabs.forEach(tab => {
    tab.addEventListener('shown.bs.tab', function(e) {
      const newPane = document.querySelector(this.getAttribute('href'));
      if (newPane) {
        const cards = newPane.querySelectorAll('.content-card');
        const links = newPane.querySelectorAll('.index-box a');
        cards.forEach(c => c.classList.remove('active'));
        links.forEach(l => l.classList.remove('active'));
        if (cards[0]) {
          cards[0].classList.add('active');
          if (links[0]) links[0].classList.add('active');
        }
      }
    });
  });

  // 🔹 2. INTERRUPTOR DE MODO OSCURO
  const themeToggle = document.getElementById('theme-toggle');
  themeToggle.addEventListener('click', () => {
    document.body.classList.toggle('dark-mode');
    if (document.body.classList.contains('dark-mode')) {
      themeToggle.innerHTML = "☀️ Modo Claro";
    } else {
      themeToggle.innerHTML = "🌓 Modo Oscuro";
    }
  });

  // 🔹 3. MANEJADOR DE DESCARGA PDF COMPLETAMENTE INTEGRAL
  const pdfButton = document.getElementById('pdf-download');
  pdfButton.addEventListener('click', () => {
    // Al hacer clic, el navegador abrirá el diálogo nativo "Guardar como PDF"
    // Gracias a las reglas '@media print' del CSS, el motor del navegador expandirá
    // automáticamente todas las pestañas y tarjetas en el documento final.
    window.print();
  });
  

});
</script>