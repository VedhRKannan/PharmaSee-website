import { Analytics } from "@vercel/analytics/react"

export const metadata = {
  title: "Pharmasee",
  description: "ADMET predictor by Pharmasee",
};

export default function RootLayout({ children }) {
  return (
    <html lang="en">
      <body>{children}</body>
    </html>
  );
}
